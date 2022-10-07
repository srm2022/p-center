package besm; //the package name and the location of the java files in your machine should match
//associated GitHub repository (at the time of writing): https://github.com/srm2022/p-center
//please amend the global variable "path" to point to the data folder on your machine

import java.io.File;
import java.io.FileNotFoundException;
import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.util.Scanner;

public class b_pnc {

    static String path = "data_for_pNextCenter\\"; //please change this to the path to folder where you put the files pmed1.txt to pmed8.txt
    static int MY_INF = 1000000;
    static int dmax;
    static int n, p, q;

    public static void main(String[] args) {
        ThreadMXBean tbbb = ManagementFactory.getThreadMXBean();
        double all_runs_f_avg, starttime, all_restarts_best_f_found_time, this_run_time, all_runs_time_avg, all_runs_time_std;
        //by t1 and f1 we mean the time and the OV obtained when reaching/exceeding the OV obtained by López Sánchez et al. (2019). 
        //by t2 and f2 we mean the time_to_best and OV obtained when run for the same amount of time as spent by López Sánchez et al. (2019).
        int this_run_f1, this_run_f2;
        double this_run_t1, this_run_t2, all_runs_t1_avg, all_runs_t1_std, all_runs_f1_avg, all_runs_f1_std, all_runs_t2_avg, all_runs_t2_std, all_runs_f2_avg, all_runs_f2_std;
        int this_restart_best_f, all_restarts_best_f;
        int all_runs_best_f; //just to record when new BKOVs are found

        int[] BKOV_targets = new int[132], OV_2019_targets = new int[132]; //these are, respectively, the BKOVs and the objective values obtained in López Sánchez et al. (2019). They only differe for two instances.
        double[] normalised_time_2019 = new double[132];
        int[] pmeds = new int[132];
        int[] ns = new int[132];
        int[] ps = new int[132];

        fillin_files_and_ns_and_ps_and_targets(BKOV_targets,OV_2019_targets, normalised_time_2019, pmeds, ns, ps);

        int i, j, t;
        long c;
        int tmpSwap;
        double tmpH;
        int tmpf;
        boolean flgImproved;

        boolean MODE_BLS = false; //true iff the experiments for the BLS in Section 4.1.1, Table 1.
        int BLS_time_out = 10 * 1000 * 1000;//100 * 1000;//millis
        int cnt_BLS_exceed_100s;

        boolean MODE_SECTION412 = true; //true iff the experiments for Section 4.1.2 in the paper
        int cnt_t1_worse_than_t_2019, cnt_f2_worse_than_f_2019;

        int run_per_instance = 10; //run per instance
        int timeout = 10 * 1000 * 000; //millisec. This is deliberately a large value so the algorithms (except BLS) run with no time limit in practice
        for (int IE = 1; IE <= 1 ; IE++) { //when IE = 1, the flat-move strategy is enabled. Use IE = 0 for BLS (A1), A3 and A5 algorithms in the paper
            System.out.println("+----+ IE = " + IE + " ----");
            for (int K = 3; K <= 3; K += 2) { //K is the parameter in the heuristic function h_K. 
                //When K is 1, the heuristic will be equal to the objective function f.
                System.out.println("-------- K = " + K + " --------");

                for (t = 0; t < 132; t++) {//for each instance
                    n = ns[t];
                    p = ps[t];
                    q = n - p;
                    int[][] D = new int[1000][1000];//[n][n];
                    int[] P = new int[p];//lists of centres
                    int[] Q = new int[q];//lists of non-centres vertices

                    //read in the instance into D
                    prepare_D_and_dmax(D, path + "pmed" + pmeds[t] + ".txt");
                    c = 2 * dmax + 1;                    

                    all_runs_f_avg = all_runs_time_avg = all_runs_time_std = 0;
                    cnt_BLS_exceed_100s = 0;
                    all_runs_t1_avg = all_runs_t1_std = all_runs_f1_avg = all_runs_f1_std = all_runs_t2_avg = all_runs_t2_std = all_runs_f2_avg = all_runs_f2_std = cnt_t1_worse_than_t_2019 = cnt_f2_worse_than_f_2019 = 0;
                            
                    all_runs_best_f = MY_INF;
                    for (int run = 0; run < run_per_instance; run++) {
                        starttime = tbbb.getCurrentThreadCpuTime() / 1000000;
                        all_restarts_best_f_found_time = tbbb.getCurrentThreadCpuTime() / 1000000 - starttime;
                        all_restarts_best_f = this_restart_best_f = this_run_f2 = dmax;
                        this_run_t1 = this_run_t2 = -1; //dummy
                        this_run_f1 = this_run_f2 = MY_INF;
                        
                        //use the following termination condition for the experiments of Table 1 in the paper. Please note that for BLS, the conditions commented out should also be used. 
                        //while (all_restarts_best_f > BKOV_targets[t] /*&& tbbb.getCurrentThreadCpuTime() / 1000000 - starttime < timeout */
                        //                                        /*&& !(IE ==0 && K == 1 && tbbb.getCurrentThreadCpuTime() / 1000000 - starttime > BLS_time_out)*/ ) {
                        //use the following termination condition for the experiments of Table 2:
                        while (tbbb.getCurrentThreadCpuTime() / 1000000 - starttime < 0.99 * normalised_time_2019[t] ||  all_restarts_best_f > OV_2019_targets[t]) {
                            initRandSol(P, Q);

                            double h = calcH(D, P, Q, K);

                            //the main loop
                            flgImproved = true;
                            while (flgImproved) {
                                flgImproved = false;
                                for (i = 0; i < p; i++) {
                                    for (j = 0; j < q; j++) {
                                        if (P[i] == Q[j]) {
                                            ee(" error ---- should not get here");
                                        }
                                        //try the neighbour
                                        tmpSwap = P[i];
                                        P[i] = Q[j];
                                        Q[j] = tmpSwap;

                                        tmpH = calcH(D, P, Q, K);
                                        if (tmpH < h) {
                                            h = tmpH;
                                            flgImproved = true;
                                            if (MODE_SECTION412) {
                                                tmpf = (int) (h / Math.pow(c, K - 1));
                                                if (tmpf <= OV_2019_targets[t] && this_run_t1 == -1) {
                                                    this_run_t1 = tbbb.getCurrentThreadCpuTime() / 1000000 - starttime;
                                                    this_run_f1 = tmpf;
                                                }
                                                if (tmpf < this_run_f2 && tbbb.getCurrentThreadCpuTime() / 1000000 - starttime < 0.99 * normalised_time_2019[t]) {
                                                    this_run_t2 = tbbb.getCurrentThreadCpuTime() / 1000000 - starttime;
                                                    this_run_f2 = tmpf;
                                                }
                                            }
                                        } else {
                                            if (tmpH > h) {
                                                tmpSwap = P[i];
                                                P[i] = Q[j];
                                                Q[j] = tmpSwap;
                                            } else {
                                                if (IE == 1) {
                                                    //no action, so the flat move is accepted
                                                } else { //undo the flat move:
                                                    tmpSwap = P[i];
                                                    P[i] = Q[j];
                                                    Q[j] = tmpSwap;
                                                }
                                            }
                                        }
                                    }//for j
                                }//for i
                            }//while flgImproved
                            this_restart_best_f = calcF(D, P, Q);
                            if (this_restart_best_f < all_restarts_best_f) {
                                all_restarts_best_f = this_restart_best_f;
                                all_restarts_best_f_found_time = tbbb.getCurrentThreadCpuTime() / 1000000 - starttime;                                
                            }
                        }//while target or timeout

                        if (MODE_BLS && tbbb.getCurrentThreadCpuTime() / 1000000 - starttime > BLS_time_out){ //this is for BLS runs for Table 1
                            this_run_time = BLS_time_out;
                            cnt_BLS_exceed_100s++;
                        }else{
                            //if (all_restarts_best_f > BKOV_targets[t]) {
                            //   System.out.println("!!! BKOV_target not met !!!");//this is just for additional information
                            //}
                            //this_run_time =  tbbb.getCurrentThreadCpuTime()/1000000  - starttime;
                            this_run_time = all_restarts_best_f_found_time;//for Table 1, it should not matter whether to use this or the above line
                        }
                        all_runs_time_avg += this_run_time;
                        all_runs_time_std += this_run_time * this_run_time;
                        all_runs_f_avg += all_restarts_best_f;    
                        if(all_restarts_best_f < all_runs_best_f){
                            all_runs_best_f = all_restarts_best_f;
                        }

                        if (MODE_SECTION412) {
                            if(this_run_t1 > normalised_time_2019[t]){
                               cnt_t1_worse_than_t_2019++;
                            }
                            if(this_run_f2 > OV_2019_targets[t]){
                               cnt_f2_worse_than_f_2019++;
                            }
                            all_runs_t1_avg += this_run_t1;
                            all_runs_t1_std += this_run_t1 * this_run_t1;
                            all_runs_f1_avg += this_run_f1;
                            all_runs_f1_std += this_run_f1 * this_run_f1;
                            all_runs_t2_avg += this_run_t2;
                            all_runs_t2_std += this_run_t2 * this_run_t2;
                            all_runs_f2_avg += this_run_f2;
                            all_runs_f2_std += this_run_f2 * this_run_f2;
                            //System.out.println(this_run_t1 + "\t" + this_run_f1 + "\t" + this_run_t2 + "\t" + this_run_f2);
                        }
                    }//for run
                    all_runs_f_avg /= run_per_instance;
                    all_runs_time_avg /= run_per_instance;
                    all_runs_time_std = Math.sqrt(all_runs_time_std / run_per_instance - all_runs_time_avg * all_runs_time_avg);
                    if (MODE_SECTION412) {
                        all_runs_t1_avg /= run_per_instance;
                        all_runs_t1_std = Math.sqrt(all_runs_t1_std / run_per_instance - all_runs_t1_avg * all_runs_t1_avg);
                        all_runs_f1_avg /= run_per_instance;
                        all_runs_f1_std = Math.sqrt(all_runs_f1_std / run_per_instance - all_runs_f1_avg * all_runs_f1_avg);

                        all_runs_t2_avg /= run_per_instance;
                        all_runs_t2_std = Math.sqrt(all_runs_t2_std / run_per_instance - all_runs_t2_avg * all_runs_t2_avg);
                        all_runs_f2_avg /= run_per_instance;
                        all_runs_f2_std = Math.sqrt(all_runs_f2_std / run_per_instance - all_runs_f2_avg * all_runs_f2_avg);
                        
                        //change from milli to seconds and round to two decimals only:
                        all_runs_t1_avg = Math.round(all_runs_t1_avg /= 10.0) / 100.0;
                        all_runs_t1_std = Math.round(all_runs_t1_std /= 10.0) / 100.0;
                        all_runs_t2_avg = Math.round(all_runs_t2_avg /= 10.0) / 100.0;
                        all_runs_t2_std = Math.round(all_runs_t2_std /= 10.0) / 100.0;
                    }
                    
                    //change from milli to seconds and round to two decimals only
                    all_runs_time_avg = Math.round(all_runs_time_avg /= 10.0) / 100.0;
                    all_runs_time_std = Math.round(all_runs_time_std /= 10.0) / 100.0;
    
                    if (MODE_BLS) {
                        System.out.println(all_runs_f_avg + "\t" + all_runs_time_avg + "\t" + all_runs_time_std + "\t" + cnt_BLS_exceed_100s);//for BLS in Table 1
                    } else {
                        if (MODE_SECTION412) {
                            //NOTE: this prints some extra results, in addition to those reported in the paper
                            System.out.println(all_runs_t1_avg + "\t" + all_runs_t1_std + "\t" + all_runs_f1_avg + "\t" + all_runs_f1_std + "\t" + (run_per_instance - cnt_t1_worse_than_t_2019) + "\t" + 
                                               all_runs_t2_avg + "\t" + all_runs_t2_std + "\t" + all_runs_f2_avg + "\t" + all_runs_f2_std + "\t" + (run_per_instance - cnt_f2_worse_than_f_2019) + "\t" + 
                                               all_runs_best_f); //for Table 2
                        } else {
                            System.out.println(all_runs_f_avg + "\t" + all_runs_time_avg + "\t" + all_runs_time_std);//for Table 1 except BLS
                        }
                    }
                    if (t == 39 || t == 73 || t == 107) {
                        System.out.println();//for readability only
                    }
                }//for each instance t
            }//for K
        }//for IE
    }//main

    //--------------
    static double calcH(int[][] D, int[] P, int[] Q, int K) //calcualtes and returns the resulting heuristic value based on the current P and U
    {
        int i, i1, i2, j, v1, v2, min, k;
        int min1, minboth, tmp1, tmpboth;

        int[] tmpd = new int[n]; //tmpd[i] shows the total cost (primary + secondary) for vi

        //first, calc distances for centres
        for (i1 = 0; i1 < p; i1++) {
            v1 = P[i1];
            min = MY_INF; //for infinity
            for (i2 = 0; i2 < p; i2++) {
                if (i2 != i1) {
                    v2 = P[i2];
                    if (D[v1][v2] < min) {
                        min = D[v1][v2];
                    }
                }
            }
            tmpd[v1] = min;
        }

        //then, calc distances for non-centres
        for (j = 0; j < q; j++) {
            v1 = Q[j];
            min1 = minboth = MY_INF;
            for (i = 0; i < p; i++) {
                v2 = P[i];
                tmp1 = D[v1][v2];
                tmpboth = tmp1 + tmpd[v2];
                if ((tmp1 < min1) || (tmp1 == min1 && tmpboth < minboth)) {
                    min1 = tmp1;
                    minboth = tmpboth;
                }
            }
            tmpd[v1] = minboth;
        }

        //calculate h based on the distances calculated above
        double[] maxH = new double[K];
        int x, k1, k2, k3;

        for (k = 0; k < K; k++) {
            maxH[k] = tmpd[0];
        }
        for (k1 = 1; k1 < n; k1++) {
            x = tmpd[k1];
            k2 = 0;
            while (k2 < K && x <= maxH[k2]) {
                k2++;
            }
            if (k2 < K) {
                for (k3 = K - 1; k3 > k2; k3--) {
                    maxH[k3] = maxH[k3 - 1];
                }
                maxH[k3] = x;
            }
        }

        double ret_val = maxH[0];
        long c = 2 * dmax + 1;
        for (k = 1; k < K; k++) {
            ret_val = ret_val * c + maxH[k];
            if (ret_val < 0) {
                ee("---error   overflow------");
            }
        }
        return ret_val;
    }//calcH
    //-----------

    static int calcF(int[][] D, int[] P, int[] Q) //calcualtes and return the resulting objective function based on the current P and U
    {
        int i, i1, i2, j, v1, v2, min, maxd;
        int min1, minboth, tmp1, tmpboth;

        maxd = 0;
        int[] tmpd = new int[n]; //tmpd[i] shows the total cost (primary + secondary) for vi

        //first, calc distances for centres
        for (i1 = 0; i1 < p; i1++) {
            v1 = P[i1];
            min = MY_INF; //for infinity
            for (i2 = 0; i2 < p; i2++) {
                if (i2 != i1) {
                    v2 = P[i2];
                    if (D[v1][v2] < min) {
                        min = D[v1][v2];
                    }
                }
            }
            tmpd[v1] = min;
            if (min > maxd) {
                maxd = min;
            }
        }

        //then, calc distances for non-centres
        for (j = 0; j < q; j++) {
            v1 = Q[j];
            min1 = minboth = MY_INF;
            for (i = 0; i < p; i++) {
                v2 = P[i];
                tmp1 = D[v1][v2];
                tmpboth = tmp1 + tmpd[v2];
                if ((tmp1 < min1) || (tmp1 == min1 && tmpboth < minboth)) {
                    min1 = tmp1;
                    minboth = tmpboth;
                }
            }
            tmpd[v1] = minboth;
            if (minboth > maxd) {
                maxd = minboth;
            }
        }

        return maxd;
    }//calcF
    //--------

    static void initRandSol(int[] P, int[] Q)//only populates P and U as a valid candidate solution
    {
        int i, j, k;
        boolean[] tmpSelected = new boolean[n];
        double tmpPbyN = p * 1.0 / n;
        for (j = 0; j < n; j++) {
            tmpSelected[j] = false;
        }

        k = 0;
        i = (int) (Math.random() * n); //random starting point
        while (k < p) {
            if (!tmpSelected[i]) {
                if (Math.random() < tmpPbyN) {
                    P[k] = i;
                    tmpSelected[i] = true;
                    k++;
                }
            }
            i = (i + 1) % n;
        }

        //copy all unselected ones into U
        k = 0;
        for (j = 0; j < n; j++) {
            if (!tmpSelected[j]) {
                Q[k] = j; //or tmpQ[j], the same
                k++;
            }
        }
    }//initRandSol
    //--------------

    static void ee(String message) {
        System.out.println(message);
        System.exit(1);
    }//ee
    //----------

    static void prepare_D_and_dmax(int[][] D, String theRawFilename) {
        int i, j, k, w;
        int tmpN = 0, tmpM = 0;  //dummy init

        //first read the number of nodes and edges:
        try {
            File infile = new File(theRawFilename);
            Scanner myReader = new Scanner(infile);
            String tmpFirstLine = myReader.nextLine().trim();
            String[] tmpParameters = new String[3];
            tmpParameters = tmpFirstLine.split(" ");
            tmpN = Integer.valueOf(tmpParameters[0]);
            tmpM = Integer.valueOf(tmpParameters[1]);
            myReader.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            ee("An error occurred.");
        }

        //Then, read the graph      
        try {
            File infile = new File(theRawFilename);
            Scanner myReader = new Scanner(infile);
            myReader.nextLine().trim();//skip the first line as the parameters have already been read

            for (i = 0; i < tmpN; i++) {
                D[i][i] = 0;
                for (j = i + 1; j < tmpN; j++) {
                    D[i][j] = D[j][i] = MY_INF; //for infinity
                }
            }

            int tmpCntM = 0;
            String e;
            String[] ijw;
            while (myReader.hasNextLine()) {
                e = myReader.nextLine().trim();
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
                ee("---- error reading the file : inconsistent number of edges----");
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            ee("An error occurred.");
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

        for (i = 0; i < n; i++) {
            D[i][i] = 0;
        }

        dmax = 0;
        for (i = 0; i < n; i++) {
            for (j = i + 1; j < n; j++) {
                if (D[i][j] > dmax) {
                    dmax = D[i][j];
                }
            }
        }
    }//prepare_D_and_dmax
    //------------

    static void fillin_files_and_ns_and_ps_and_targets(int[] BKOV_targets, int[] OV_2019_targets, double[] normalised_time_2019, int[] pmeds, int[] ns, int[] ps) {
        int t;
        String tmps;
        String[] tmpsplits;

        //pmedd
        for (t = 0; t < 10; t++) {
            pmeds[t] = 1;
        }
        for (t = 10; t < 20; t++) {
            pmeds[t] = 2;
        }
        for (t = 20; t < 30; t++) {
            pmeds[t] = 3;
        }
        for (t = 30; t < 40; t++) {
            pmeds[t] = 4;
        }

        for (t = 40; t < 57; t++) {
            pmeds[t] = 1;
        }
        for (t = 57; t < 74; t++) {
            pmeds[t] = 2;
        }

        for (t = 74; t < 91; t++) {
            pmeds[t] = 3;
        }
        for (t = 91; t < 108; t++) {
            pmeds[t] = 4;
        }

        for (t = 108; t < 116; t++) {
            pmeds[t] = 6;
        }
        for (t = 116; t < 124; t++) {
            pmeds[t] = 7;
        }
        for (t = 124; t < 132; t++) {
            pmeds[t] = 8;
        }

        //ns
        tmps = "10,20,20,30,30,40,40,40,50,50,10,20,20,30,30,40,40,40,50,50,10,20,20,30,30,40,40,40,50,50,10,20,20,30,30,40,40,40,50,50,"
                + "60,60,60,70,70,70,80,80,80,90,90,90,90,100,100,100,100,60,60,60,70,70,70,80,80,80,90,90,90,90,100,100,100,100,"
                + "60,60,60,70,70,70,80,80,80,90,90,90,90,100,100,100,100,60,60,60,70,70,70,80,80,80,90,90,90,90,100,100,100,100,"
                + "150,150,150,150,200,200,200,200,150,150,150,150,200,200,200,200,150,150,150,150,200,200,200,200";
        tmpsplits = tmps.split(",");
        for (t = 0; t < 132; t++) {
            ns[t] = Integer.valueOf(tmpsplits[t]);
        }

        //ps
        tmps = "5,5,10,5,10,5,10,20,10,20,5,5,10,5,10,5,10,20,10,20,5,5,10,5,10,5,10,20,10,20,5,5,10,5,10,5,10,20,10,20,"
                + "10,20,30,10,20,30,10,20,30,10,20,30,50,10,20,30,50,10,20,30,10,20,30,10,20,30,10,20,30,50,10,20,30,50,"
                + "10,20,30,10,20,30,10,20,30,10,20,30,50,10,20,30,50,10,20,30,10,20,30,10,20,30,10,20,30,50,10,20,30,50,"
                + "20,30,50,80,20,30,50,80,20,30,50,80,20,30,50,80,20,30,50,80,20,30,50,80";
        tmpsplits = tmps.split(",");
        for (t = 0; t < 132; t++) {
            ps[t] = Integer.valueOf(tmpsplits[t]);
        }
        
        //BKOV_targets
        tmps = "84,120,95,126,95,144,111,89,110,89,121,147,99,169,110,164,112,96,140,99,77,145,77,157,122,157,105,77,125,87,126,139,125,173,122,175,122,85,126,91," +
               "112,91,89,119,99,73,133,105,91,133,108,91,70,133,108,97,74,140,99,96,138,102,96,138,109,97,140,109,97,96,135,109,96,96," +
               "124,97,73,121,97,82,121,93,86,148,105,93,93,151,113,93,93,135,93,79,146,102,85,146,114,91,147,112,92,82,147,119,96,82," +
               "79,71,62,56,79,72,68,54,69,62,59,59,73,68,63,52,74,61,58,58,84,77,68,68";
        tmpsplits = tmps.split(",");
        for (t = 0; t < 132; t++) {
            BKOV_targets[t] = Integer.valueOf(tmpsplits[t]);
            OV_2019_targets[t] = BKOV_targets[t];
        }
        OV_2019_targets[8] = 111;
        OV_2019_targets[10] = 128;
        
        //normalised_time_2019
        tmps = "0.02,0.04,0.12,0.09,0.28,0.24,0.57,1.67,0.94,3.22,0.01,0.05,0.1,0.1,0.36,0.2,0.56,1.68,0.84,3.28,0.01,0.06,0.1,0.1,0.32,0.2,0.56,1.7,0.86,2.77,0.01,0.04,0.11,0.09,0.29,0.19,0.54,1.73,0.93,3.41," +
               "1.26,5.08,11.5,1.84,7.74,18.85,2.44,11.64,29.16,3.34,15.05,38.29,103.71,4.36,19.32,51.51,143.32,1.43,6.1,13.61,1.97,8.83,20.65,2.78,12.4,30.84,3.56,16.52,41.52,104.57,3.97,18.81,50.26,145.13," +
               "1.47,5.14,10.83,2.15,8.66,19.58,2.89,12.04,29.93,2.82,12.6,33.52,89.92,3.56,17.13,44.94,124.77,1.46,5.43,13.89,2.04,8.43,19.17,2.81,12.05,29.39,3.63,16.47,43.52,102.93,4.71,21.9,57.82,154.18," +
               "33.86,77.19,200.17,479.01,49.78,150.29,493.55,1151.12,22.9,66.1,206.33,415.03,44.05,129.25,505.03,1213.66,23.76,64.84,194.59,430.29,41.21,122.14,394.95,1140.05";
        tmpsplits = tmps.split(",");
        for (t = 0; t < 132; t++) {
            normalised_time_2019[t] = 1000 * Double.valueOf(tmpsplits[t]) / 2.60; //2.60 is the CPU performance ratio based on CPU Benchmark (3027 vs. 1166). 1000 is multiplied for millisecs
                                                                           
        }
    }//fillin_files_and_ns_and_ps_and_targets
    //------------

}//end class
