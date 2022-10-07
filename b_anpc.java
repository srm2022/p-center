package besm; //the package name and the location of the java files in your machine should match
//associated GitHub repository (at the time of writing): https://github.com/srm2022/p-center
//please amend the global variable "path" to point to the  data folder on your machine

import java.io.File;
import java.io.FileNotFoundException;
import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.util.Scanner;

public class b_anpc {
    static String path = "pmed_data\\"; //please change this to the path to pmed data on your machine

    static int a = 2;
    static int MY_INF = 1000000;

    static String[] files;
    static int[] ns, ps;
    static int[] targets;

    static double dmax;
    static double c;
    static int n, p;
    static double[][] D;
    static int[] C; 
    static int nC;
    static double f;

    static int[][] N; 
    static int[][] Ninverse; 
    static boolean[] x;

    static long[][] tb, tbt;
    static int[] tb_index_u, tb_index_p;
    static int ntb;
    static long all_runs_move_cnt;

    public static void main(String[] args) {
        ThreadMXBean tbbb = ManagementFactory.getThreadMXBean();
        fillin_files_and_ns_and_ps_and_targets();

        int t;
        double starttime, this_run_time, all_runs_time_avg, all_runs_time_std;

        int[][] F; //keeps list of facilities assigned to each client

        int v, i, j, pi, uj, y, tmpcntcurpos, tmpcntp, inxi, fa, tmpp, tmpq, curpos, cc;
        double h, lsbesth, tmph, leastworseh, this_restart_best_f, all_restarts_best_f;
        boolean[] this_restart_best_x, all_restarts_best_x;
        boolean flgmoved;

        int nL;

        int[] tmpP, tmpQ, tmpuinverse;

        int MAX_ITR = 10; 
        int timeout = 10000000; 
        for (int IE = 0; IE <= 1; IE++) {//IE = 1 means IE strategy is enabled.
            for (int K = 1; K <= 3; K += 2) {//the parameter of h_K
                System.out.println("------ IE = " + IE + " ------ K = " + K);
                for (t = 0; t < 40; t++) {
                    n = ns[t];
                    p = ps[t];
                    prepare_D_and_dmax(t);
                    c = 1 * dmax + 1; //or any number greater than this
                    tmpP = new int[n];
                    tmpQ = new int[n];  
                    tmpuinverse = new int[n];
                    C = new int[n];
                    nC = 0;

                    N = new int[n][n];
                    Ninverse = new int[n][n];
                    fillin_s_and_Ninverse();

                    x = new boolean[n];
                    F = new int[n][a + 1];

                    tb = new long[n][n];
                    tbt = new long[n][n];
                    tb_index_u = new int[n * n];
                    tb_index_p = new int[n * n];
                    ntb = 0;
                    all_runs_move_cnt = 0;

                    for (i = 0; i < n; i++) {
                        for (j = 0; j < n; j++) {
                            tb[i][j] = tbt[i][j] = 0;
                        }
                    }

                    //Li and Luj implement SCL in the pseudocode in the paper
                    int[] Li = new int[(n - p) * p];
                    int[] Luj = new int[(n - p) * p];

                    this_restart_best_x = new boolean[n];
                    all_restarts_best_x = new boolean[n];

                    all_runs_time_avg = all_runs_time_std = 0;
                    for (int run = 0; run < MAX_ITR; run++) {

                        starttime = tbbb.getCurrentThreadCpuTime()/1000000;
                        all_restarts_best_f = dmax;
                        while (Math.round(all_restarts_best_f * 100.0) / 100 > targets[t] && tbbb.getCurrentThreadCpuTime()/1000000 - starttime < timeout) {
                            h = initRandSol(x, F, K);
                            for (v = 0, tmpp = tmpq = 0; v < n; v++) {
                                if (x[v]) {
                                    tmpP[tmpp] = v;
                                    tmpuinverse[v] = -1;
                                    tmpp++;
                                } else {
                                    tmpQ[tmpq] = v;
                                    tmpuinverse[v] = tmpq;
                                    tmpq++;
                                }
                            }
                            this_restart_best_f = f;
                            for (int ix = 0; ix < n; this_restart_best_x[ix] = x[ix], ix++); //copy x to this_restart_best_x
                            if (this_restart_best_f < all_restarts_best_f) {
                                all_restarts_best_f = this_restart_best_f;
                                for (int ix = 0; ix < n; all_restarts_best_x[ix] = this_restart_best_x[ix], ix++);
                            }

                            //LS
                            reset_tabu();

                            lsbesth = h;
                            flgmoved = true;
                            long lsstarttime = tbbb.getCurrentThreadCpuTime()/1000000;
                            while (flgmoved && Math.round(this_restart_best_f * 100.0) / 100 > targets[t] && tbbb.getCurrentThreadCpuTime()/1000000 - lsstarttime < timeout) {
                                flgmoved = false;
                                nL = 0;
                                leastworseh = Double.MAX_VALUE;
                                cc = C[(int) (Math.random() * nC)];
                                curpos = Ninverse[cc][F[cc][a - 1]];
                                for (y = (int) (Math.random() * curpos), tmpcntcurpos = 0; !flgmoved && tmpcntcurpos < curpos; y = (y + 1) % curpos, tmpcntcurpos++) {
                                    uj = N[cc][y];
                                    if (x[uj]) {
                                        continue;
                                    }
                                    if (D[cc][uj] < f) {
                                        //find a p 
                                        for (inxi = (int) (Math.random() * p), tmpcntp = 0; !flgmoved && tmpcntp < p; inxi = (inxi + 1) % p, tmpcntp++) {
                                            fa = tmpP[inxi];
                                            tmph = h_value_if_swapped(uj, fa, x, F, K);
                                            if (tmph < h || tmph == h && IE > 0 && (tbt[uj][fa] + tb[uj][fa] <= all_runs_move_cnt)) {
                                                flgmoved = true;
                                                all_runs_move_cnt++;
                                                pi = tmpP[inxi];

                                                h = tmpswap(uj, pi, inxi, x, F, K, tmpQ, tmpP, tmpuinverse);
                                                if (h < lsbesth) {
                                                    lsbesth = h;
                                                    reset_tabu();
                                                } else {
                                                    entabu(uj, pi);
                                                }

                                                if (f < this_restart_best_f) {
                                                    this_restart_best_f = f;
                                                    for (int ix = 0; ix < n; this_restart_best_x[ix] = x[ix], ix++);
                                                }
                                            }//if better than oe equal to f
                                            else {
                                                if (tbt[uj][fa] + tb[uj][fa] <= all_runs_move_cnt) {
                                                    if (tmph < leastworseh) {
                                                        leastworseh = tmph;
                                                        Li[0] = inxi;
                                                        Luj[0] = uj;
                                                        nL = 1;
                                                    //else { //this part and the related few lines of code are commented out because SCL is not needed really due to the random order of selecting (uj,inxi) pairs
                                                    //    if (tmph == leastworseh) {
                                                    //        Li[nL] = inxi;
                                                    //        Luj[nL] = uj;
                                                    //        nL++;
                                                    //    }
                                                    }
                                                }
                                            }//if-else not moved
                                        }//for p
                                    }//if beneficial
                                }//for y

                                if (!flgmoved) {
                                    if (nL > 0) {
                                        i = (int) (Math.random() * nL);
                                        inxi = Li[i];
                                        uj = Luj[i];
                                        pi = tmpP[inxi];
                                        h = tmpswap(uj, pi, inxi, x, F, K, tmpQ, tmpP, tmpuinverse);
                                        flgmoved = true;
                                        all_runs_move_cnt++;
                                        entabu(uj, pi);
                                    }else{
                                        //System.out.println("------------ nL = 0 ---------------");
                                    }//nL
                                }//if should do the best possile move
                            }//LS

                            if (this_restart_best_f < all_restarts_best_f) {
                                all_restarts_best_f = this_restart_best_f;
                                for (int ix = 0; ix < n; all_restarts_best_x[ix] = this_restart_best_x[ix], ix++);
                            }
                        }//outer while not target

                        if (all_restarts_best_f != calcF(all_restarts_best_x)) {
                            ee("final answer is wrong!");
                        }

                        this_run_time = (tbbb.getCurrentThreadCpuTime()/1000000 - starttime);
                        all_runs_time_avg += this_run_time;
                        all_runs_time_std += this_run_time * this_run_time;
                    }// for run
                    all_runs_time_avg /= MAX_ITR;
                    all_runs_time_std = Math.sqrt(all_runs_time_std / MAX_ITR - all_runs_time_avg * all_runs_time_avg);
                    //change from milli to seconds and rounds to two decimals only
                    all_runs_time_avg = Math.round(all_runs_time_avg /= 10.0)/100.0;
                    all_runs_time_std = Math.round(all_runs_time_std /= 10.0)/100.0;
                    System.out.println(all_runs_time_avg+ "\t" + all_runs_time_std);
                    if (t == 39 || t == 83 || t == 123) {
                        System.out.println();//for readability only
                    }
                }//for each instance t
            }//IE          
        }//mode
    }//main

    //--------------
    static double tmpswap(int uj, int pi, int inxi, boolean[] x, int[][] F, int K, int[] tmpU, int[] tmpP, int[] tmpuinverse) {
        double fnew;

        fnew = swap(uj, pi, x, F, K);

        int inxj = tmpuinverse[uj];
        tmpP[inxi] = uj;
        tmpuinverse[uj] = -1;
        tmpU[inxj] = pi;
        tmpuinverse[pi] = inxj;

        return fnew;
    }//tmpswap
    //---------

    static void entabu(int uj, int pi) {
        if (tb[uj][pi] == 0)//not already in the indexes
        {
            tb_index_u[ntb] = uj;
            tb_index_p[ntb] = pi;
            ntb++;
        }
        tb[uj][pi] = tb[pi][uj] = Long.MAX_VALUE / 1000000;
    }//entabu
    //---------

    static void reset_tabu() {
        for (int k = 0; k < ntb; k++) {
            tb[tb_index_u[k]][tb_index_p[k]] = tb[tb_index_p[k]][tb_index_u[k]] = 0;
            tbt[tb_index_u[k]][tb_index_p[k]] = tbt[tb_index_p[k]][tb_index_u[k]] = 0;
        }
        ntb = 0;
    }//reset_tabu
    //---------

    static double h_value_if_swapped(int come, int gone, boolean[] x, int[][] F, int K) {

        double tmpf = -1;
        double[] tmpd = new double[n];
        int gonepos, comepos, lastpos;

        for (int v = 0; v < n; v++) {
            if (!x[v] && v != come || v == gone) {
                gonepos = Ninverse[v][gone];
                comepos = Ninverse[v][come];
                lastpos = Ninverse[v][F[v][a - 1]];
                if (gonepos > lastpos) {
                    if (comepos > lastpos) {
                        tmpd[v] = D[v][F[v][a - 1]];
                    } else {
                        if (a == 1) {
                            tmpd[v] = D[v][come];
                        } else {
                            tmpd[v] = Math.max(D[v][F[v][a - 2]], D[v][come]);
                        }
                    }
                } else {
                    if (comepos > Ninverse[v][F[v][a]]) {
                        tmpd[v] = D[v][F[v][a]];
                    } else {
                        if (a == 1) {
                            tmpd[v] = D[v][come];
                        } else {
                            if (gone == F[v][a - 1]) {
                                tmpd[v] = Math.max(D[v][F[v][a - 2]], D[v][come]);
                            } else {
                                tmpd[v] = Math.max(D[v][F[v][a - 1]], D[v][come]);
                            }
                        }
                    }
                }

                if (tmpd[v] > tmpf) {
                    tmpf = tmpd[v];
                }
            }
        }
        if (K == 1) {
            return tmpf;
        } else {
            x[gone] = false;
            x[come] = true;
            double tmph = h(tmpd, K);
            x[gone] = true;
            x[come] = false;
            return tmph;
        }
    }//f_value_if_swapped
    //---------

    static double h(double[] tmpd, int K) {
        //calc h based on the d values
        double[] maxH = new double[K];
        int k, k2, k3;
        double z;

        for (k = 0; k < K; k++) {
            maxH[k] = 0;
        }
        for (int v = 0; v < n; v++) {//starting from 0 not 1
            if (!x[v]) {
                z = tmpd[v];
                k2 = 0;
                while (k2 < K && z <= maxH[k2]) {
                    k2++;
                }
                if (k2 < K) {
                    for (k3 = K - 1; k3 > k2; k3--) {
                        maxH[k3] = maxH[k3 - 1];
                    }
                    maxH[k3] = z;
                }
            }
        }

        return calc_h_for_list(maxH, K);
    }//h
    //-----------

    static double calc_h_for_list(double[] list, int nhist) {
        double ret_val = list[0];

        for (int k = 1; k < nhist; k++) {
            ret_val = ret_val * c + list[k];
        }
        return ret_val;
    }//calc_h_for_list
    //-----------

    static double swap(int come, int gone, boolean[] x, int[][] F, int K) {
        int k, k2;
        double[] tmpd = new double[n];
        nC = -1;
        f = -1;
        x[come] = true;
        x[gone] = false;

        for (int v = 0; v < n; v++) {
            //add
            k = 0;
            while (k <= a && Ninverse[v][come] > Ninverse[v][F[v][k]]) {
                k++;
            }
            if (k <= a) {
                for (k2 = a; k2 > k; k2--) {
                    F[v][k2] = F[v][k2 - 1];
                }
                F[v][k2] = come;
            }

            //remove
            k = 0;
            while (k <= a && gone != F[v][k]) {
                k++;
            }
            if (k <= a) {
                for (k2 = k; k2 < a; k2++) {
                    F[v][k2] = F[v][k2 + 1];
                }
                F[v][a] = N[v][nextpos(v, Ninverse[v][F[v][a]])];
            }

            if (!x[v]) {
                tmpd[v] = D[v][F[v][a - 1]];
                if (tmpd[v] > f) {
                    f = tmpd[v];
                    C[0] = v;
                    nC = 1;
                } else {
                    if (tmpd[v] == f) {
                        C[nC] = v;
                        nC++;
                    }
                }
            }
        }
        if (K == 1) {
            return f;
        } else {
            return h(tmpd, K);
        }
    }//swap
    //--------

    static double initRandSol(boolean[] x, int[][] F, int K) {
        int k, fa, v, minpos;
        double mind;
        double[] tmpd = new double[n];

        double tmpPbyN = p * 1.0 / n;
        for (fa = 0; fa < n; fa++) {
            x[fa] = false;
        }

        k = 0;
        fa = (int) (Math.random() * n); //random starting point
        while (k < p) {
            if (!x[fa]) {
                if (Math.random() < tmpPbyN) {
                    x[fa] = true;
                    k++;
                }
            }
            fa = (fa + 1) % n;
        }

        nC = -1;
        f = -1;

        for (v = 0; v < n; v++) {
            minpos = -1;
            for (k = 0; k < a + 1; k++) {
                minpos = nextpos(v, minpos);
                F[v][k] = N[v][minpos];
            }//for k

            if (!x[v]) {
                mind = D[v][F[v][a - 1]];
                if (mind > f) {
                    f = mind;
                    C[0] = v;
                    nC = 1;
                } else {
                    if (mind == f) {
                        C[nC] = v;
                        nC++;
                    }
                }

                tmpd[v] = mind;
            }
        }
        if (K == 1) {
            return f;
        } else {
            return h(tmpd, K);
        }
    }//initRandSol
    //----------

    static int nextpos(int c, int curpos) {
        int pos = 1 + curpos;
        while (!x[N[c][pos]]) {
            pos++;
        }
        return pos;
    }//nextpos
    //----------

    static void ee(String message) {
        System.out.println(message);
        System.exit(1);
    }//ee
    //----------

    static double calcF(boolean[] x)//returns the objective value
    {
        int v1, v2;
        double tmp, maxd;
        double[] mins = new double[a];

        maxd = -1;

        //calc ds for non-centres
        for (v1 = 0; v1 < n; v1++) {
            if (x[v1] == false) {// simply remove this requirement for the sister problem, the a-all neighbour problem defined in Khuller et al. (2000)- no further changes is needed
                for (int k = 0; k < a; k++) {
                    mins[k] = MY_INF;
                }
                for (v2 = 0; v2 < n; v2++) {
                    if (x[v2]) {
                        tmp = D[v1][v2];
                        addtolist(mins, a, tmp, true);
                    }
                }
                if (mins[a - 1] > maxd) {
                    maxd = mins[a - 1];
                }
            }
        }
        return maxd;
    }//calcF
    //--------

    static void addtolist(double[] list, int num, double element, boolean asc) {
        int k, k2;

        if (asc) {
            k = 0;
            for (k = 0; k < num && element >= list[k]; k++);
            if (k < num) {
                for (k2 = num - 2; k2 >= k; list[k2 + 1] = list[k2], k2--);
                list[k] = element;
            }
        } else {//i.e. descendingly
            k = 0;
            for (k = 0; k < num && element <= list[k]; k++);
            if (k < num) {
                for (k2 = num - 2; k2 >= k; list[k2 + 1] = list[k2], k2--);
                list[k] = element;
            }
        }
    }//addtolist
    //--------

    static void fillin_s_and_Ninverse() {
        int j, c;

        for (c = 0; c < n; c++) {
            for (j = 0; j < n; j++) {
                N[c][j] = j;
            }
            quicksort(c, 0, n - 1);
        }

        //now Ninverse
        for (c = 0; c < n; c++) {
            for (j = 0; j < n; j++) {
                Ninverse[c][N[c][j]] = j;
            }
        }
    }//fillin_s_and_Ninverse
    //----------

    static void quicksort(int c, int from, int to) {
        int i, j, tmp;
        double pivotitem;

        if (from < to) {
            pivotitem = D[N[c][from]][c];
            j = from;
            for (i = from + 1; i <= to; i++) {
                if (D[N[c][i]][c] < pivotitem) {
                    j++;
                    tmp = N[c][i];
                    N[c][i] = N[c][j];
                    N[c][j] = tmp;
                }
            }

            tmp = N[c][j];
            N[c][j] = N[c][from];
            N[c][from] = tmp;

            quicksort(c, from, j - 1);
            quicksort(c, j + 1, to);
        }
    }//quicksort
    //----------

    static boolean prepare_D_and_dmax(int t) {
        if (t < 40) {
            String theRawFilename;

            theRawFilename = path + files[t];

            int i, j, k, w, tmp_inf;
            int tmpN = 0, tmpM = 0, tmpP = 0;//dummy init
            //read in the graph      
            try {
                File infile = new File(theRawFilename);
                Scanner myReader = new Scanner(infile);
                String tmpFirstLine = myReader.nextLine().trim();
                String[] tmpParameters = new String[3];
                tmpParameters = tmpFirstLine.split(" ");
                tmpN = Integer.valueOf(tmpParameters[0]);
                tmpM = Integer.valueOf(tmpParameters[1]);
                tmpP = Integer.valueOf(tmpParameters[2]);

                tmp_inf = Integer.MAX_VALUE / 2; //for infinity; devided by 2 to avoid overflow when adding two of them in Floyd alg.
                D = new double[tmpN][tmpN];
                for (i = 0; i < tmpN; i++) {
                    D[i][i] = 0;
                    for (j = i + 1; j < tmpN; j++) {
                        D[i][j] = D[j][i] = tmp_inf;
                    }
                }

                int tmpCntM = 0;
                String e;
                String[] ijw;
                while (myReader.hasNextLine()) {
                    e = myReader.nextLine().trim();
                    //System.out.println(e);
                    ijw = e.split(" ");
                    i = Integer.valueOf(ijw[0]);
                    j = Integer.valueOf(ijw[1]);
                    w = Integer.valueOf(ijw[2]);
                    if ((i == j) || (w <= 0) || (i > tmpN) || (j > tmpN)) {
                        ee("error: inconsistency in the input file!");
                    } else {
                        D[i - 1][j - 1] = D[j - 1][i - 1] = w;
                        tmpCntM++;
                    }
                }
                myReader.close();
                if (tmpM != tmpCntM) {
                    ee("error: inconsistent m");
                }
            } catch (FileNotFoundException e) {
                System.out.println("error reading input file!");
                e.printStackTrace();
                System.exit(1);
            }

            //run Floyd
            for (k = 0; k < tmpN; k++) {
                for (i = 0; i < tmpN; i++) {
                    for (j = 0; j < tmpN; j++) {
                        if (D[i][k] + D[k][j] < D[i][j]) {
                            D[i][j] = D[i][k] + D[k][j];
                        }
                    }
                }
            }

            n = tmpN;
            p = tmpP;

            dmax = 0;
            for (i = 0; i < n; i++) {
                for (j = i + 1; j < n; j++) {
                    if (D[i][j] > dmax) {
                        dmax = D[i][j];
                    }
                }
            }

            if (n != ns[t] || p != ps[t]) {
                ee("inconsistency in readng data files!");
            }
        } else {
            String path2 = "tsp_data\\"; 
            String theRawFilename;

            theRawFilename = path2 + files[t];

            int i, j, k, tmpN;
            Double x, y, dx, dy, d;
            Double[] tmpx, tmpy;
            //read in the graph      
            try {
                File infile = new File(theRawFilename);
                Scanner myReader = new Scanner(infile);
                //skip first 3 lines:
                myReader.nextLine();
                myReader.nextLine();
                myReader.nextLine();
                String tmpDimensionLine = myReader.nextLine().trim();
                String[] tmps = new String[3];
                tmps = tmpDimensionLine.split(" ");
                if (tmps[0].charAt(tmps[0].length() - 1) == ':') {
                    tmpN = Integer.valueOf(tmps[1]);
                } else {
                    tmpN = Integer.valueOf(tmps[2]);
                }
                if (tmpN != ns[t]) {
                    myReader.close();
                    ee("inconsistent number of nodes in the files!");
                }
                //skip hte othe two lines:
                myReader.nextLine();
                myReader.nextLine();

                tmpx = new Double[tmpN];
                tmpy = new Double[tmpN];

                int tmpcntrow = 0;
                String e;
                String[] kxy;
                e = myReader.nextLine().trim();
                while (e.charAt(0) != '1') {
                    e = myReader.nextLine().trim();
                }

                while (e.compareTo("EOF") != 0) {
                    tmpcntrow++;
                    kxy = e.split(" ");
                    k = Integer.valueOf(kxy[0]);
                    x = Double.valueOf(kxy[1]);
                    y = Double.valueOf(kxy[2]);
                    if ((k != tmpcntrow) || (tmpcntrow > tmpN)) {
                        ee("error: inconsistency in the input file!");
                    } else {
                        tmpx[k - 1] = x;
                        tmpy[k - 1] = y;
                    }
                    e = myReader.nextLine().trim();
                }
                myReader.close();

                //fill in D
                D = new double[tmpN][tmpN];
                dmax = -1;
                for (i = 0; i < tmpN; i++) {
                    D[i][i] = 0;
                    for (j = i + 1; j < tmpN; j++) {
                        dx = tmpx[i] - tmpx[j];
                        dy = tmpy[i] - tmpy[j];
                        d = Math.sqrt(dx * dx + dy * dy);
                        D[i][j] = D[j][i] = Math.round(d * 100.0) / 100.0;
                        if (d > dmax) {
                            dmax = d;
                        }
                    }
                }
            } catch (FileNotFoundException e) {
                System.out.println("error reading input file");
                e.printStackTrace();
                System.exit(1);
            }
        }//if else t
        return true;
    }//prepare_D_and_dmax
    //----------

    static void fillin_files_and_ns_and_ps_and_targets() {
        files = new String[166];
        ns = new int[166];
        ps = new int[166];
        targets = new int[166];

        int k;
        String tmps;
        String[] tmpsplits;

        //pmeds
        for (k = 0; k < 40; k++) {
            files[k] = "pmed" + (k + 1) + ".txt";
        }
        for (k = 0; k < 34; k++) {
            ns[k] = (k / 5 + 1) * 100;
        }
        ns[34] = ns[35] = ns[36] = 800;
        ns[37] = ns[38] = ns[39] = 900;

        tmps = "5,10,10,20,33,5,10,20,40,67,5,10,30,60,100,5,10,40,80,133,5,10,50,100,167,5,10,60,120,200,5,10,70,140,5,10,80,5,10,90";
        tmpsplits = new String[40];
        tmpsplits = tmps.split(",");
        for (k = 0; k < 40; k++) {
            ps[k] = Integer.valueOf(tmpsplits[k]);
        }

        //small TSPs
        for (k = 40; k < 44; k++) {
            files[k] = "pr226" + ".tsp";
            ns[k] = 226;
        }
        for (k = 44; k < 48; k++) {
            files[k] = "pr264" + ".tsp";
            ns[k] = 264;
        }
        for (k = 48; k < 52; k++) {
            files[k] = "pr299" + ".tsp";
            ns[k] = 299;
        }
        for (k = 52; k < 56; k++) {
            files[k] = "pr439" + ".tsp";
            ns[k] = 439;
        }
        for (k = 56; k < 60; k++) {
            files[k] = "pcb442" + ".tsp";
            ns[k] = 442;
        }
        for (k = 60; k < 64; k++) {
            files[k] = "kroA200" + ".tsp";
            ns[k] = 200;
        }
        for (k = 64; k < 68; k++) {
            files[k] = "kroB200" + ".tsp";
            ns[k] = 200;
        }
        for (k = 68; k < 72; k++) {
            files[k] = "lin318" + ".tsp";
            ns[k] = 318;
        }
        for (k = 72; k < 76; k++) {
            files[k] = "gr202" + ".tsp";
            ns[k] = 202;
        }
        for (k = 76; k < 80; k++) {
            files[k] = "d493" + ".tsp";
            ns[k] = 493;
        }
        for (k = 80; k < 84; k++) {
            files[k] = "d657" + ".tsp";
            ns[k] = 657;
            ps[k] = (k - 83) * 10;
        }
        int tmpp = 40;
        for (k = 40; k < 84; k++) {
            ps[k] = tmpp;
            if (tmpp == 5) {
                tmpp = 40;
            } else {
                tmpp /= 2;
            }
        }

        //large TSPs
        for (k = 84; k < 99; k++) {
            files[k] = "u1060" + ".tsp";
            ns[k] = 1060;
            ps[k] = (k - 83) * 10;
        }
        for (k = 99; k < 109; k++) {
            files[k] = "rl1323" + ".tsp";
            ns[k] = 1323;
            ps[k] = (k - 98) * 10;
        }
        for (k = 109; k < 124; k++) {
            files[k] = "u1817" + ".tsp";
            ns[k] = 1817;
            ps[k] = (k - 108) * 10;
        }

        //for that pcb files used in memtics2008 but not in grasp2017:
        for (k = 124; k < 138; k++) {
            files[k] = "pcb3038.tsp";
            ns[k] = 3038;
            if (k < 129) {
                ps[k] = (k - 123) * 10;
            } else {
                ps[k] = (k - 128) * 50;
            }
        }

        //other TSPs for a-neighbor in the future esA
        for (k = 138; k < 142; k++) {
            files[k] = "att48.tsp";
            ns[k] = 48;
            ps[k] = (k - 137) * 10;
        }
        for (k = 142; k < 152; k++) {
            files[k] = "eil101.tsp";
            ns[k] = 101;
            ps[k] = (k - 141) * 10;
        }
        for (k = 152; k < 166; k++) {
            files[k] = "ch150.tsp";
            ns[k] = 150;
            ps[k] = (k - 151) * 10;
        }

        tmps = "150,121,121,97,63,99,80,70,49,28,68,60,43,34,23,52,45,34,25,19,45,44,27,20,15,43,36,22,17,13,34,33,19,14,34,31,19,33,26,16";
        tmpsplits = tmps.split(",");
        for (int t = 0; t < 40; t++) {
            targets[t] = (int) (Math.round(Double.valueOf(tmpsplits[t]) * 100.0) / 100);
        }
    }//fillin_files_and_ns_and_ps_and_targets
    //----------
}//end class


