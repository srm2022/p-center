package besm; //the package name and the location of the java files in your machine should match
//associated GitHub repository (at the time of writing): https://github.com/srm2022/p-center
//please amend the global variable "path" to point to the data folder on your machine

//Before using this source code, please change the global variables "path1" and "path2" to point to the right data folders on your machine.
import java.io.File;
import java.io.FileNotFoundException;
import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.util.Arrays;
import java.util.Scanner;

public class b_pc {

    static String path1 = "pmed_data//"; //please change this to the path to pmed data folder on your machine.
    static String path2 = "tsp_data//"; //and this to your TSP data folder.

    static int MY_INF = 1000000;

    static String[] files;
    static int[] ns, ps;
    static double[] targets, targets_pbs, targets_grasppr;

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
    static int[] P, Q, Qinverse;
    static double[] m_p;
    static int[][] m_p_v;
    static int[] m_p_v_cnt;

    static long[][] tb, tbt;
    static int[] tb_index_u, tb_index_p;
    static int ntb;
    static long[] tabuP, tabuPt;
    static int[] tabuP_index;
    static int ntabuP;
    static long all_runs_move_cnt;

    static int[] F;
    static int[] Finx; //to indicate the index of each client in  Z
    static int[] F2;
    static int[] F2inx; //to indicate the index of each client in  Z2
    static int[][] Z;
    static int[][] Z2;
    static int[] nZ;
    static int[] nZ2;
    
    static double epsilon = 0.001; //NOTE: it assumes maximum two decimal figures in the distance values; otherwise, please amend this value.

    static int IE, UNFLAT;
    static Boolean DEBUG_MODE = false;

    public static void main(String[] args) {
  
        
   ThreadMXBean tbbb = ManagementFactory.getThreadMXBean();
   
        fillin_files_and_ns_and_ps_and_targets();

        int t;
        double all_runs_f_avg, starttime, all_restarts_best_f_found_time, this_run_time, all_runs_time_avg, all_runs_time_std, target_pbs_f, target_grasppr_f, target_pbs_t, target_grasppr_t;

        int k, i, j, pi, uj, y, tmpcntcurpos, tmpcntp, inxi, fa, curpos, cc, dfu, nC1;
        double r, h, lsbesth, minm, this_restart_best_f, all_restarts_best_f, tmp_value, tmp1 = 0, tmp2 = 0, tmp_swap_time = 0, tmp_other_time = 0;
        double[] promising_outputs = new double[3];
        boolean[] all_restarts_best_x;
        boolean[] this_restart_best_x;
        boolean flgmoved;

        //These implement SCL in the pseudocode in the paper:
        int BEST_VAL_CNT = 3;
        double[] best_value_list = new double[BEST_VAL_CNT]; 
        int[] best_value_uj = new int[BEST_VAL_CNT];
        int[] best_value_inxi = new int[BEST_VAL_CNT]; 

        boolean TARGET_PBS_REPORT = false;
        boolean TARGET_GRASPPR_REPORT = false;
        int run_per_instance = 100;//10; //run per instance
        int timeout = 10 * 1000 * 1000; //millisec.
        for (IE = 1; IE <= 1; IE++) { //IE = 1 means IE strategy is enabled. Use 0 for algorithms without it.
            for (UNFLAT = 1; UNFLAT <= 1; UNFLAT++) { //unflat = 1 means unflattening is enabled. Use 0 for algorithms without it.
                System.out.println("------ IE = " + IE + " ------------ unflat = " + UNFLAT);
                double all_instances_all_runs_time = tbbb.getCurrentThreadCpuTime()/1000000;
                for (t = 0; t < 40; t++) {
                    if(DEBUG_MODE)
                        System.out.print("t = " + t + " : -------------------------------------------");
                    n = ns[t];
                    p = ps[t];
                    prepare_D_and_dmax(t);
                    P = new int[n];
                    Q = new int[n];
                    Qinverse = new int[n];
                    C = new int[n];
                    nC = 0;
                    m_p = new double[n];
                    m_p_v = new int[n][n];
                    m_p_v_cnt = new int[n];

                    N = new int[n][n];
                    Ninverse = new int[n][n];
                    fillin_N_and_Ninverse();

                    x = new boolean[n];
                    F = new int[n];
                    Finx = new int[n];
                    F2 = new int[n];
                    F2inx = new int[n];
                    Z = new int[n][n];
                    Z2 = new int[n][n];
                    nZ = new int[n];
                    nZ2 = new int[n];

                    tb = new long[n][n];
                    tbt = new long[n][n];
                    tb_index_u = new int[n * n];
                    tb_index_p = new int[n * n];
                    ntb = 0;
                    all_runs_move_cnt = 0;
                    tabuP = new long[n];
                    tabuPt = new long[n];
                    tabuP_index = new int[n];
                    ntabuP = 0;

                    for (i = 0; i < n; i++) { 
                        for (j = 0; j < n; j++) {
                            tb[i][j] = tbt[i][j] = -1; //not 0;
                        }
                    }
                    for (i = 0; i < n; i++) {
                        tabuP[i] = tabuPt[i] = -1; //not 0;
                    }

                    all_restarts_best_x = new boolean[n];
                    this_restart_best_x = new boolean[n];

                    all_runs_f_avg = all_runs_time_avg = all_runs_time_std = 0;
                    for (int run = 0; run < run_per_instance; run++) {
                        if(DEBUG_MODE)
                            System.out.print("\n - run = " + run);
                        starttime =  tbbb.getCurrentThreadCpuTime()/1000000;
                        all_restarts_best_f_found_time =  tbbb.getCurrentThreadCpuTime()/1000000;
                        all_restarts_best_f = dmax;
                        this_restart_best_f = -1;//dummy

                        target_pbs_f = target_grasppr_f = target_pbs_t = target_grasppr_t = -1;

                        double ls_while_itr_coeff = 1;
                        boolean first_time = true;
                        while (Math.round(all_restarts_best_f * 100.0) / 100.0 > targets[t] &&  tbbb.getCurrentThreadCpuTime()/1000000 - starttime < timeout) {
                            if(DEBUG_MODE)
                              System.out.print(" - init ");
                            if (first_time || Math.random() < 0.33) {
                                h = reset_solution(UNFLAT, null, null);
                                first_time = false;
                            }else{
                                if(Math.random() < 0.5){
                                    h = reset_solution(UNFLAT, this_restart_best_x, null);
                                }else{
                                    h = reset_solution(UNFLAT, null, all_restarts_best_x);
                                }                                
                            }

                            lsbesth = h;
                            if (DEBUG_MODE) {
                                test_data(x);
                                tmp_test_m_p_stuff_and_nZ();
                            }
                            this_restart_best_f = f;
                            for (int ix = 0; ix < n; this_restart_best_x[ix] = x[ix], ix++);
                            if (f < all_restarts_best_f) {//though unlikely
                                all_restarts_best_f = f;
                                all_restarts_best_f_found_time =  tbbb.getCurrentThreadCpuTime() / 1000000 - starttime;
                                for (int ix = 0; ix < n; all_restarts_best_x[ix] = x[ix], ix++); //this line copies array x to all_restarts_best_x.
                            }

                            //LS
                            reset_tabu();

                            flgmoved = true;
                            int ls_while_itr = 0;
                            while (flgmoved && ls_while_itr < (int)(ls_while_itr_coeff * (n-p) * p) && Math.round(this_restart_best_f * 100.0) / 100.0 > targets[t] &&  tbbb.getCurrentThreadCpuTime()/1000000 - starttime < timeout) {
                                ls_while_itr++;
                                flgmoved = false;
                                minm = Double.MAX_VALUE;
                                for (k = 0; k < BEST_VAL_CNT; k++) {
                                    best_value_list[k] = Double.MAX_VALUE;
                                    best_value_uj[k] = -1;
                                    best_value_inxi[k] = -1;
                                }
                                
                                cc = C[(int) (Math.random() * nC)];
                                curpos = Ninverse[cc][F[cc]];

                                for (y = (int) (Math.random() * curpos), tmpcntcurpos = 0; !flgmoved && tmpcntcurpos < curpos; y = (y + 1) % curpos, tmpcntcurpos++) {
                                    uj = N[cc][y];
                                    if (D[cc][uj] < f) {
                                        //find a p 
                                        for (inxi = (int) (Math.random() * p), tmpcntp = 0; !flgmoved && tmpcntp < p; inxi = (inxi + 1) % p, tmpcntp++) {
                                            fa = P[inxi];

                                            if(promising2(uj, fa, minm, promising_outputs))//NOTE: function promising1 may be used instead (which is simpler, with not too different performance)
                                            {
                                                dfu = (int) promising_outputs[0];
                                                nC1 = (int) promising_outputs[1]; //NOTE: nC1 does not have a valid value unless dfu = 0 and must not be used otherwise
                                                tmp_value = promising_outputs[2];
                                                if(dfu >= 0 && tmp_value < 0 || tmp_value > minm)
                                                    ee("inconsistency here");

                                                if ((dfu < 0 || dfu == 0 && nC1 < nC && UNFLAT == 1)
                                                             || dfu == 0 && (UNFLAT == 0 || nC1 == nC) && IE > 0) {
                                                    flgmoved = true;
                                                    all_runs_move_cnt++;
                                                    pi = P[inxi];

                                                    tmp1 = tbbb.getCurrentThreadCpuTime() / 1000;
                                                    tmp_other_time += (tmp1 - tmp2);
                                                    h = tmpswap(uj, pi, inxi);
                                                    tmp2 = tbbb.getCurrentThreadCpuTime() / 1000;
                                                    tmp_swap_time += (tmp2 - tmp1);

                                                    if (dfu == 0 && (UNFLAT == 0 || nC1 == nC)) {//IE must have been 1 and no need to explicitly include
                                                        entabu(uj, pi);
                                                    }
                                                    if (h < lsbesth) {
                                                        lsbesth = h;
                                                        ls_while_itr = 0;
                                                        reset_tabu();
                                                    }
                                                    if (f < this_restart_best_f) {
                                                        this_restart_best_f = f;
                                                        for (int ix = 0; ix < n; this_restart_best_x[ix] = x[ix], ix++);
                                                        if (f < all_restarts_best_f) {
                                                            all_restarts_best_f = f;
                                                            all_restarts_best_f_found_time = tbbb.getCurrentThreadCpuTime() / 1000000 - starttime;
                                                            for (int ix = 0; ix < n; all_restarts_best_x[ix] = x[ix], ix++);
                                                            //System.out.println("all_restarts_best_f = " + all_restarts_best_f + "( " + Math.round(all_restarts_best_f_found_time/10.0)/100.0 + " )");
                                                            if (target_pbs_t == -1 && targets_pbs[t] != targets[t] && f <= targets_pbs[t]) {
                                                                target_pbs_f = f;
                                                                target_pbs_t = tbbb.getCurrentThreadCpuTime() / 1000000 - starttime;
                                                            }
                                                            if (target_grasppr_t == -1 && targets_grasppr[t] != targets[t] && f <= targets_grasppr[t]) {
                                                                target_grasppr_f = f;
                                                                target_grasppr_t = tbbb.getCurrentThreadCpuTime() / 1000000 - starttime;
                                                            }
                                                        }
                                                    }
                                                }//if better than or equal to f
                                                else {
                                                    //update SCL
                                                    for (k = 0; k < BEST_VAL_CNT && tmp_value >= best_value_list[k]; k++)
                                                        ;
                                                    if (k < BEST_VAL_CNT) {
                                                        for (int k2 = BEST_VAL_CNT - 1; k2 > k; k2--) {
                                                            best_value_list[k2] = best_value_list[k2 - 1];
                                                            best_value_uj[k2] = best_value_uj[k2 - 1];
                                                            best_value_inxi[k2] = best_value_inxi[k2 - 1];
                                                        }
                                                        best_value_list[k] = tmp_value;
                                                        best_value_uj[k] = uj;
                                                        best_value_inxi[k] = inxi;
                                                    }
                                                    minm = best_value_list[BEST_VAL_CNT - 1];
                                                }//if-else not moved
                                            }
                                        }//for p
                                    }//if beneficial
                                }//for y

                                if (!flgmoved) {
                                    if (best_value_inxi[0] >= 0) {
                                        for (k = 0; k < BEST_VAL_CNT && best_value_inxi[k] >= 0; k++)
                                            ;
                                        if(k == 0)
                                            ee("invalid best_value_inxi: it should have had at least one valid member!");
                                        r = Math.random() * k * (k + 1) / 2;
                                        int k2 = 0;
                                        r -= k;
                                        while (r > 0) {
                                            k2++;
                                            r -= (k - k2);
                                        }
                                        inxi = best_value_inxi[k2];
                                        uj = best_value_uj[k2];
                                        pi = P[inxi];

                                        tmp1 = tbbb.getCurrentThreadCpuTime() / 1000;
                                        tmp_other_time += (tmp1 - tmp2);
                                        h = tmpswap(uj, pi, inxi);
                                        tmp2 = tbbb.getCurrentThreadCpuTime() / 1000;
                                        tmp_swap_time += (tmp2 - tmp1);

                                        flgmoved = true;
                                        all_runs_move_cnt++;
                                        entabu(uj, pi);
                                        entabuP(uj);
                                    }
                                }
                                if (!flgmoved) {
                                    if(DEBUG_MODE)
                                        System.out.println(" SCL empty! ");
                                }
                            }//LS
                            ls_while_itr_coeff *= 1.1;
                        }//outer while not target

                        if (all_restarts_best_f != calcF(all_restarts_best_x)) {
                            ee("final solution incorrect!");
                        }
                        this_run_time = ( tbbb.getCurrentThreadCpuTime()/1000000 - starttime);
                        all_runs_time_avg += this_run_time;
                        all_runs_time_std += this_run_time * this_run_time;

                        if(DEBUG_MODE){
                        //if(true){
                            System.out.print(all_restarts_best_f + "\t" + this_run_time / 1000 + "\t");
                            if(TARGET_PBS_REPORT && targets_pbs[t] != targets[t]){
                                target_pbs_t = Math.round(target_pbs_t /= 10.0) / 100.0;
                                System.out.print("(" + target_pbs_f + "@" + target_pbs_t + ")\t");
                            }
                            if(TARGET_GRASPPR_REPORT && targets_grasppr[t] != targets[t]){
                                target_grasppr_t = Math.round(target_grasppr_t /= 10.0) / 100.0;
                                System.out.print("(" + target_grasppr_f + "@" + target_grasppr_t + ")\t");
                            }
                        }
                        all_runs_f_avg += all_restarts_best_f;
                    }// for run
                    all_runs_f_avg /= run_per_instance;
                    all_runs_time_avg /= run_per_instance;
                    all_runs_time_std = Math.sqrt(all_runs_time_std / run_per_instance - all_runs_time_avg * all_runs_time_avg);
                    //change from milli to seconds and round to two decimals only
                    all_runs_time_avg = Math.round(all_runs_time_avg /= 10.0) / 100.0;
                    all_runs_time_std = Math.round(all_runs_time_std /= 10.0) / 100.0;
                    System.out.println(all_runs_f_avg + "\t" + all_runs_time_avg + "\t" + all_runs_time_std);
                    if (t == 39 || t == 83 || t == 123) {
                        System.out.println();//for readability only
                    }
                }//for each instance t
                all_instances_all_runs_time =  tbbb.getCurrentThreadCpuTime()/1000000 - all_instances_all_runs_time;
                System.out.println("all_instances_all_runs_time = " + all_instances_all_runs_time / 1000);
            }//
        }//  
        System.out.println("total tmp_other_time in millis = " + tmp_other_time/1000 + ", total tmp_swap_time in millis= " + tmp_swap_time/1000);

    }//main
    //-----------------------

    static boolean promising1(int uj, int fa, double minm, double[] promising_outputs) {
        //This function returns false if the potential swap (uj, fa) is either tabued or is no better than those in SCL.
        //If it returns true, it also prepared output values in promising_outputs array to be used by the main function.
        //This function is much simpler than promising2.

        int k, tmpv, tmpcc, dfu, nC1, tmpcnt, tmpF2;
        double tmpd, m, tmp_value;
        double[] m_value_outputs = new double[2];

        if(minm < f)
            ee("inconsistency: minm should not be less than f");
        
        if ( (tbt[uj][fa] + tb[uj][fa] > all_runs_move_cnt) || (tabuPt[fa] + tabuP[fa] > all_runs_move_cnt) ){
            return false;
        }

        tmp_value = -1;//-1 means not needed, i.e. when dfu = -1 or this function returns false
        nC1 = -1; //-1 means not needed, i.e. when dfu != 0

        int count1 = 0;
        for (k = 0; k < nC; k++) {
            tmpcc = C[k];
            tmpd = D[tmpcc][uj];
            if (F[tmpcc] != fa && tmpd >= f) {
                count1++;
            }
        }

        int nA = 0;            
        m = 0;
        tmpcnt = nZ[fa];
        for (k = 0; k < tmpcnt && m < minm; k++) {
            tmpv = Z[fa][k];
            tmpF2 = F2[tmpv];
            tmpd = D[tmpv][uj];
            if (D[tmpv][tmpF2] < tmpd) {
                tmpd = D[tmpv][tmpF2];
            }
            if (tmpd > m) {
                m = tmpd;
            }
            if (tmpd == f) {
                nA++;
            }
        }

        if (m >= minm) {
            return false;
        }

        if (m < f && count1 == 0) {
            dfu = -1;
        } else {
            if (m == f || m < f && count1 > 0) {
                dfu = 0;
                nC1 = nA + count1;
                tmp_value = f;
            } else {//i.e. m > f
                dfu = 1;
                tmp_value = m;
            }
        }
        
        if (DEBUG_MODE) {
            //double-check with an alternative (previous) logic:
            int dfu2;
            m_value_if_swapped(fa, uj, m_value_outputs);
            double m2 = m_value_outputs[0];
            int nC12 = 0;
            if (m2 <= f) {
                if (m2 == f) {
                    nC12 = (int) m_value_outputs[1];
                }
                if (nC > 1) {
                    for (int k2 = 0; k2 < nC; k2++) {
                        tmpcc = C[k2];
                        tmpd = D[tmpcc][uj];
                        if (F[tmpcc] != fa && tmpd >= f) {
                            nC12++;
                        }
                    }
                }
                if (nC12 > 0) {
                    dfu2 = 0;
                } else {
                    dfu2 = -1;
                }
            } else {
                dfu2 = 1;
            }

            if (dfu != 10 && (dfu != dfu2 || dfu == 0 && dfu2 == 0 && nC12 != nC1)) {
                ee("inconsistency 2");
            }
        }
        promising_outputs[0] = dfu;
        promising_outputs[1] = nC1;
        promising_outputs[2] = tmp_value;

        return true;
    }//promising1
    //----------

    static boolean promising2(int uj, int fa, double minm, double[] promising_outputs) {
        //This function may be ignored and the function promising1 may be used instead. This function may be slightly faster in some cases but no significant difference was observed. 
        //It uses some temporary variables to avoid the cost of calculating m in some cases. These variables are:
        //count1: the number of critical vertices that are not assigned to fa and whose costs (f) is no more than their distance to uj (i.e. they remain critical even if uj becomes a facility).
        //m_LB: a lower bound(LB) on the actual value of m. For example, if this value is no less than minm, then m cannot be so either, hence we will not need to put the input par (uj, fa) in SCL.
        //calculation of m_LB is based on a list (m_p_v[fa]) of nodes assigned to fa whose F2 values (kept in m_p[fa]) are the maximum (among those assigned to fa). The number of such nodes is kept in m_p_v_cnt[fa]. These variables are updated once a move is performed. m_LB is the maximum of the costs of the node in m_p_v[fa] if the facility fa is replaced with uj. 
        //If you use promising1 instead of this function, please remove m_p_v[], m_p_v_cnt[] and m_p[] from the program.
        
        int k, tmpv, tmpcc, dfu, nC1, tmpcnt, tmpF2;
        double tmp_m_p, tmpd, tmpdouble, m, tmp_value;
        double[] m_value_outputs = new double[2];

        if(minm < f)
            ee("inconsistency: minm should not be less than f");
        
        if ( (tbt[uj][fa] + tb[uj][fa] > all_runs_move_cnt) || (tabuPt[fa] + tabuP[fa] > all_runs_move_cnt) ){
            return false;
        }

        tmp_value = -1;//-1 means not needed, i.e. when dfu = -1 or this function returns false
        nC1 = -1; //-1 means not needed, i.e. when dfu != 0


        int count1 = 0;//count 1 is the number of critical clients not assigned to fa such that their costs does not decrease by making uj a new facility
        for (k = 0; k < nC; k++) {
            tmpcc = C[k];
            tmpd = D[tmpcc][uj];
            if (F[tmpcc] != fa && tmpd >= f) {
                count1++;
            }
        }

        //cases:
        tmp_m_p = m_p[fa];
        if (tmp_m_p < f) {
            if (count1 == 0) {
                dfu = -1;
            } else {
                dfu = 0;
                nC1 = count1;
                tmp_value = f;
            }
        } else {
            //in this case, first calculate m_LB
            double m_LB = -1;//m_LB os a lower bound (LB) on the actual value of m. 
            int cnt_m_LB = -1;//dummy init
            for (k = 0; k < m_p_v_cnt[fa]; k++) {
                tmpv = m_p_v[fa][k];
                tmpdouble = Math.min(D[tmpv][uj], tmp_m_p);
                if (tmpdouble > m_LB) {
                    m_LB = tmpdouble;
                    cnt_m_LB = 1;
                } else {
                    if (tmpdouble == m_LB) {
                        cnt_m_LB++;
                    }
                }
            }
            //check if should skip because it falls into some cases where it will be tabued 
            if (m_LB >= minm) {
                return false;
            }
            if (tmp_m_p == f) {
                if (m_LB < f && count1 == 0) {
                    dfu = -1;
                } else { //i.e. max_... ==f || cnt_... >0
                    dfu = 0;
                    nC1 = count1;
                    if (m_LB == f) {
                        nC1 += cnt_m_LB;
                    }
                    tmp_value = f;
                }
            } else {//i.e. tmp_m_p > f
                //this case may in turn be split into further subcases (but can be ignored for simplicity)
                int nA = 0;
                m = 0;
                tmpcnt = nZ[fa];
                for (k = 0; k < tmpcnt && m < minm; k++) {
                    tmpv = Z[fa][k];
                    tmpF2 = F2[tmpv];
                    tmpd = D[tmpv][uj];
                    if (D[tmpv][tmpF2] < tmpd) {
                        tmpd = D[tmpv][tmpF2];
                    }
                    if (tmpd > m) {
                        m = tmpd;
                    }
                    if (tmpd == f) {
                        nA++;
                    }
                }

                if (m >= minm) {
                    return false;
                }

                if (m < f && count1 == 0) {
                    dfu = -1;
                } else {
                    if (m == f || m < f && count1 > 0) {
                        dfu = 0;
                        nC1 = nA + count1;
                        tmp_value = f;
                    } else {//i.e. m > f
                        dfu = 1;
                        tmp_value = m;
                    }
                }
            }
        }

        if (DEBUG_MODE) {
            //double-check with an alternative (previous) logic:
            int dfu2;
            m_value_if_swapped(fa, uj, m_value_outputs);
            double m2 = m_value_outputs[0];
            int nC12 = 0;
            if (m2 <= f) {
                if (m2 == f) {
                    nC12 = (int) m_value_outputs[1];
                }
                if (nC > 1) {
                    for (int k2 = 0; k2 < nC; k2++) {
                        tmpcc = C[k2];
                        tmpd = D[tmpcc][uj];
                        if (F[tmpcc] != fa && tmpd >= f) {
                            nC12++;
                        }
                    }
                }
                if (nC12 > 0) {
                    dfu2 = 0;
                } else {
                    dfu2 = -1;
                }
            } else {
                dfu2 = 1;
            }

            if (dfu != 10 && (dfu != dfu2 || dfu == 0 && dfu2 == 0 && nC12 != nC1)) {
                ee("inconsistency 2");
            }
        }
        promising_outputs[0] = dfu;
        promising_outputs[1] = nC1;
        promising_outputs[2] = tmp_value;

        return true;
    }//promising2
    //----------
    
    static void tmp_test_m_p_stuff_and_nZ() {//for debug purpose only
        //test m_p stuff:
        for (int node = 0; node < n; node++) {
            if (x[node]) {
                int fa = node;
                double tmp_m_p = -1;
                int tmp_m_p_cnt = 0;
                int tmpnZ = 0;
                for (int v = 0; v < n; v++) {
                    if (F[v] == fa) {
                        tmpnZ++;
                        if (D[v][F2[v]] > tmp_m_p) {
                            tmp_m_p = D[v][F2[v]];
                            tmp_m_p_cnt = 1;
                        } else {
                            if (D[v][F2[v]] == tmp_m_p) {
                                tmp_m_p_cnt++;
                            }
                        }
                    }
                }
                if (tmp_m_p != m_p[fa] || tmp_m_p_cnt != m_p_v_cnt[fa] || tmpnZ != nZ[fa]) {
                    ee("inconsistency 3");
                }
            } else {
                int u = node;
                if (m_p[u] != -1 || m_p_v_cnt[u] > 0) {
                    ee("inconsistency 4");
                }
            }
        }
    }//tmp_test_m_p_stuff_and_nZ

    //---------
    static double tmpswap(int uj, int pi, int inxi) {
        double fnew;

        add_facility_and_update(uj);
        if (DEBUG_MODE) {
            tmp_test_m_p_stuff_and_nZ();
        }
        fnew = remove_facility_and_update(pi);
        if (DEBUG_MODE) {
            tmp_test_m_p_stuff_and_nZ();
        }

        int inxj = Qinverse[uj];
        P[inxi] = uj;
        Qinverse[uj] = -1;
        Q[inxj] = pi;
        Qinverse[pi] = inxj;

        if (DEBUG_MODE) {
            test_data(x);
        }

        return fnew;
    }//tmpswap
    //---------

    static void entabuP(int new_pi) {
        if (tabuP[new_pi] == -1)//not already in tabuP_index
        {
            tabuP_index[ntabuP] = new_pi;
            ntabuP++;
        }
        tabuP[new_pi] = 1;
        tabuPt[new_pi] = all_runs_move_cnt;
    }//entabuP
    //---------

    static void entabu(int uj, int pi) {
        int cap2 = (int) (0.1 * (n - p) * p);

        if (tb[uj][pi] == -1)//not already in the index
        {
            tb_index_u[ntb] = uj;
            tb_index_p[ntb] = pi;
            ntb++;
        }
        tb[uj][pi] = (tb[uj][pi] == -1) ? 1 : (int) Math.ceil(2 * tb[uj][pi]);

        if (tb[uj][pi] > cap2) {
            tb[uj][pi] = cap2;
        }
        tb[pi][uj] = tb[uj][pi];
        tbt[pi][uj] = tbt[uj][pi] = all_runs_move_cnt;
    }//entabu
    //---------

    static void reset_tabu() {
        int tmpk, pi, uj;

        for (tmpk = 0; tmpk < ntb; tmpk++) {
            uj = tb_index_u[tmpk];
            pi = tb_index_p[tmpk];
            tb[uj][pi] = tb[pi][uj] = -1; //not 0;
            tbt[uj][pi] = tbt[pi][uj] = -1; //not 0;
        }
        ntb = 0;

        for (tmpk = 0; tmpk < ntabuP; tmpk++) {
            pi = tabuP_index[tmpk];
            tabuP[pi] = tabuPt[pi] = -1; //not 0;
        }
        ntabuP = 0;
    }//reset_tabu
    //---------

    static void m_value_if_swapped(int gone, int come, double[] outputs) {
        int k, v, nA;
        double m;
        double tmpd;

        m = nA = 0;
        for (k = 0; k < nZ[gone]; k++) {
            v = Z[gone][k];
            tmpd = Math.min(D[v][come], D[v][F2[v]]);
            if (tmpd > m) {
                m = tmpd;
            }
            if (tmpd == f) {
                nA++;
            }
        }
        outputs[0] = m;
        outputs[1] = nA;
    }//m_value_if_swapped
    //---------

    static double add_facility_and_update(int come) {
        int tmpOldF, tmpOldF2;
        int fa, uj;
        double[] tmpd = new double[n];
        nC = -1;
        f = -1;
        x[come] = true;

        m_p[come] = -1;
        m_p_v_cnt[come] = 0;
        int[] tmplist = new int[n];
        int ntmplist = 0;//for tmplist which is to be the list of facilities whose m_p stuff needs revision

        for (int v = 0; v < n; v++) { //for each node
            if (Ninverse[v][come] < Ninverse[v][F[v]]) {
                tmpOldF = F[v];
                tmpOldF2 = F2[v];

                deassignF(v, tmpOldF);
                assignF(v, come);

                disconnectF2(v, tmpOldF2);
                connectF2(v, tmpOldF);

                //update p_v stuff
                //first the old one
                if (Math.abs(D[v][tmpOldF2] - m_p[tmpOldF]) < epsilon) {
                    tmplist[ntmplist] = tmpOldF;
                    ntmplist++;
                }
                //then the new one
                if (D[tmpOldF][v] > m_p[come]) {
                    m_p[come] = D[tmpOldF][v];
                    m_p_v[come][0] = v;
                    m_p_v_cnt[come] = 1;
                } else {
                    if (Math.abs(D[tmpOldF][v] - m_p[come]) < epsilon) {
                        m_p_v[come][m_p_v_cnt[come]] = v;
                        m_p_v_cnt[come]++;
                    }
                }
            } else {
                if (Ninverse[v][come] < Ninverse[v][F2[v]]) {
                    tmpOldF2 = F2[v];
                    disconnectF2(v, tmpOldF2);
                    connectF2(v, come);

                    //also m_p stuff
                    if (Math.abs(D[v][tmpOldF2] - m_p[F[v]]) < epsilon && D[v][come] <= D[v][tmpOldF2]) {
                        tmplist[ntmplist] = F[v];
                        ntmplist++;
                    }
                }
            }

            tmpd[v] = D[v][F[v]];
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

        //revise m_p for those needed
        for (int k = 0; k < ntmplist; k++) {
            fa = tmplist[k];
            m_p[fa] = -1;
            m_p_v_cnt[fa] = 0;
            for (int k2 = 0; k2 < nZ[fa]; k2++) {
                uj = Z[fa][k2];
                if (D[uj][F2[uj]] > m_p[fa]) {
                    m_p[fa] = D[uj][F2[uj]];
                    m_p_v[fa][0] = uj;
                    m_p_v_cnt[fa] = 1;
                } else {
                    if (D[uj][F2[uj]] == m_p[fa]) {
                        m_p_v[fa][m_p_v_cnt[fa]] = uj;
                        m_p_v_cnt[fa]++;
                    }
                }
            }
        }//for k

        if (UNFLAT == 0) {
            return f;
        } else {
            return n * f + nC;
        }
    }//add_facility_and_update
    //--------

    static double remove_facility_and_update(int gone) {
        int k, tmpcnt, v, tmpOldF2, tmpNewF2;
        double tmpd, tmpoldd, oldf;
        boolean flg_f_increased;
        oldf = f; //NOTE: f or nC must NOT be initialised (to, say, 0 or -1) because of the incremental nature of the function.
        x[gone] = false;//crucial to do before the actual process
        m_p[gone] = -1; //to be neat otherwise not needed
        m_p_v_cnt[gone] = 0;
        tmpcnt = nZ[gone];
        flg_f_increased = false;

        for (k = 0; k < tmpcnt; k++) {
            v = Z[gone][k];
            //NOTE: v is not de-assigned from gone at this stage; it is done at the end when all other nodes are processed
            tmpOldF2 = F2[v];
            assignF(v, tmpOldF2);

            tmpoldd = D[v][gone];
            tmpd = D[v][F[v]];
            if (tmpd > f) {
                f = tmpd;
                C[0] = v;
                nC = 1;
                flg_f_increased = true;
            } else {
                if (tmpd == f) {
                    if (tmpd < tmpoldd) {
                        ee("inconsistency 5");
                    } else {
                        if (tmpoldd > oldf) {
                            ee("inconsistency 6");
                        } else { //i.e. tmpoldd == oldf
                            if (tmpoldd < oldf || flg_f_increased) { //i.e tmpd == f && v is not in the critical list)
                                C[nC] = v;
                                nC++;
                            } else {//i.e. tmpd == f && tmpoldd == oldf
                                ;
                            }
                        }
                    }
                }
            }

            disconnectF2(v, tmpOldF2);
            tmpNewF2 = N[v][nextpos(v, Ninverse[v][F2[v]])];
            connectF2(v, tmpNewF2);

            //update m_p stuff
            if (D[v][tmpNewF2] > m_p[tmpOldF2]) {
                m_p[tmpOldF2] = D[v][tmpNewF2];
                m_p_v[tmpOldF2][0] = v;
                m_p_v_cnt[tmpOldF2] = 1;
            } else {
                if (D[v][tmpNewF2] == m_p[tmpOldF2]) {
                    m_p_v[tmpOldF2][m_p_v_cnt[tmpOldF2]] = v;
                    m_p_v_cnt[tmpOldF2]++;
                }
            }
        }
        nZ[gone] = 0; //NOTE: CRUCIAL to deassign all from gone here and NOT in the above code (as noted above)

        tmpcnt = nZ2[gone];
        for (k = 0; k < tmpcnt; k++) {
            v = Z2[gone][k];
            //disconnectF2(v, gone); //CRUCIAL to comment this out because it would otherwise change nZ[gone] used in the loops condition
            tmpNewF2 = N[v][nextpos(v, Ninverse[v][gone])];
            connectF2(v, tmpNewF2);

            //update m_p stuff
            if (D[v][tmpNewF2] > m_p[F[v]]) {
                m_p[F[v]] = D[v][tmpNewF2];
                m_p_v[F[v]][0] = v;
                m_p_v_cnt[F[v]] = 1;
            } else {
                if (D[v][tmpNewF2] < D[v][gone]) {
                    ee("inconsistency 7");
                }
                if (D[v][tmpNewF2] == m_p[F[v]] && D[v][tmpNewF2] > D[v][gone]) {
                    int tmpF = F[v];
                    int tmpmovcnt = m_p_v_cnt[F[v]];
                    if (tmpF > n - 1 || tmpF < 0 || tmpmovcnt < 0 || tmpmovcnt > n - 1) {
                        ee("error: tmpF=" + tmpF + ", tmpmovcnt=" + tmpmovcnt);
                    }
                    m_p_v[F[v]][m_p_v_cnt[F[v]]] = v;
                    m_p_v_cnt[F[v]]++;
                }
            }
        }
        nZ2[gone] = 0; //NOTE: CRUCIAL to disconnect all F2s from gone here and NOT in the above code (as commented out)

        if (DEBUG_MODE) {
            //test nC:
            int tmpnc = 0;
            for (v = 0; v < n; v++) {
                if (D[v][F[v]] == f) {
                    tmpnc++;
                }
            }
            if (tmpnc != nC) {
                ee("error: tmpnc != nC");
            }
        }

        if (UNFLAT == 0) {
            return f;
        } else {
            return n * f + nC;
        }
    }//remove_facility_and_update
    //--------

    static void assignF(int v, int newF) {
        F[v] = newF;
        Z[newF][nZ[newF]] = v;
        Finx[v] = nZ[newF];
        nZ[newF]++;
    }//connectF2
    //--------

    static void deassignF(int v, int gone) {
        int v3rd = Z[gone][nZ[gone] - 1];
        Finx[v3rd] = Finx[v];
        Z[gone][Finx[v]] = v3rd;
        nZ[gone]--;
    }//connectF2
    //--------

    static void connectF2(int v, int newF2) {
        F2[v] = newF2;
        Z2[newF2][nZ2[newF2]] = v;
        F2inx[v] = nZ2[newF2];
        nZ2[newF2]++;
    }//connectF2
    //--------

    static void disconnectF2(int v, int oldF2) {
        int v3rd = Z2[oldF2][nZ2[oldF2] - 1];
        F2inx[v3rd] = F2inx[v];
        Z2[oldF2][F2inx[v]] = v3rd;
        nZ2[oldF2]--;
    }//connectF2
    //--------

    static double reset_solution(int unflat, boolean[] this_restart_best_x, boolean [] all_restarts_best_x) {
        int k, fa, i, v, minpos, minpos2, minfa, minfa2;
        double mind, mind2;
        
        if (this_restart_best_x != null) {
            k = 0;
            for (v = 0; v < n; v++) {
                x[v] = false;
            }
            while (k < p) {
                v = (int) (Math.random() * n);
                if (x[v]) {
                    continue;
                }
                if (this_restart_best_x[v] == true) {
                    x[v] = true;
                    k++;
                }
            }
        } else {
            if (all_restarts_best_x != null) {
                k = 0;
                for (v = 0; v < n; v++) {
                    x[v] = false;
                }
                while (k < p) {
                    v = (int) (Math.random() * n);
                    if (x[v]) {
                        continue;
                    }
                    if (all_restarts_best_x[v] == true) {
                        x[v] = true;
                        k++;
                    }
                }
            } else {//then fully random restart
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
            }
        }
       
        //fill in data structures
        int tmpp, tmpq;
        for (v = 0, tmpp = tmpq = 0; v < n; v++) {
            if (x[v]) {
                P[tmpp] = v;
                Qinverse[v] = -1;
                tmpp++;
            } else {
                Q[tmpq] = v;
                Qinverse[v] = tmpq;
                tmpq++;
            }
        }
        if(tmpp != p || tmpq != n - tmpp)
            ee("error: incorrect tmpp or tmpq!");

        nC = -1;
        f = -1;

        for (i = 0; i < n; i++) { //n not p (i.e., for all potential fa)
            nZ[i] = nZ2[i] = 0;
            m_p[i] = -1;
            m_p_v_cnt[i] = 0;
        }
        for (v = 0; v < n; v++) {
            minpos = nextpos(v, -1);
            minfa = N[v][minpos];
            mind = D[minfa][v];
            minpos2 = nextpos(v, minpos);
            minfa2 = N[v][minpos2];
            mind2 = D[minfa2][v];

            F[v] = minfa;
            Z[minfa][nZ[minfa]] = v;
            Finx[v] = nZ[minfa];//excellent
            nZ[minfa]++;

            F2[v] = minfa2;
            Z2[minfa2][nZ2[minfa2]] = v;
            F2inx[v] = nZ2[minfa2];//excellent
            nZ2[minfa2]++;

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

            if (x[minfa]) {
                if (mind2 > m_p[minfa]) {
                    m_p[minfa] = mind2;
                    m_p_v[minfa][0] = v;
                    m_p_v_cnt[minfa] = 1;
                } else {
                    if (mind2 == m_p[minfa]) {
                        m_p_v[minfa][m_p_v_cnt[minfa]] = v;
                        m_p_v_cnt[minfa]++;
                    }
                }
            }
        }

        if (unflat == 0) {
            return f;
        } else {
            return n * f + nC;
        }
    }//reset_solution
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

    static double calcF(boolean[] x){//returns the objective value
        int v1, v2, cnttrue;
        double tmp, min, maxd;

        maxd = -1;
        cnttrue = 0;
        //calc ds for non-centres
        for (v1 = 0; v1 < n; v1++) {
            if (x[v1] == false) {
                min = MY_INF;
                for (v2 = 0; v2 < n; v2++) {
                    if (x[v2]) {
                        tmp = D[v1][v2];
                        if (tmp < min) {
                            min = tmp;
                        }
                    }
                }
                if (min > maxd) {
                    maxd = min;
                }
            } else {
                cnttrue++;
            }
        }
        if (cnttrue != p) {
            ee("incorrect number of open facilities!");
        }
        return maxd;
    }//calcF
    //--------

    static void fillin_N_and_Ninverse() {
        int j, v;

        for (v = 0; v < n; v++) {
            for (j = 0; j < n; j++) {
                N[v][j] = j;
            }
            quicksort(v, 0, n - 1);
        }

        //now Ninverse
        for (v = 0; v < n; v++) {
            for (j = 0; j < n; j++) {
                Ninverse[v][N[v][j]] = j;
            }
        }
    }//fillin_N_and_Ninverse
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

    static void prepare_D_and_dmax(int t) {
        if (t < 40) {
            String theRawFilename;

            theRawFilename = path1 + files[t];

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
                    ee("error: inconsistent m!");
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
                        if (D[i][j] > dmax) {
                            dmax = d;
                        }
                    }
                }
            } catch (FileNotFoundException e) {
                System.out.println("error reading input file");
                System.exit(1);
            }
        }//if else t
    }//prepare_D_and_dmax
    //----------

    static void fillin_files_and_ns_and_ps_and_targets() {
        files = new String[124];
        ns = new int[124];
        ps = new int[124];
        targets = new double[124];
        targets_pbs = new double[124];
        targets_grasppr = new double[124];
        
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

        //targets:
        tmps = "127,98,93,74,48,84,64,55,37,20,59,51,36,26,18,47,39,28,18,13,40,38,22,15,11,38,32,18,13,9,30,29,15,11,30,27,15,29,23,13,"
                + "650.00, 1365.65, 2326.48, 3720.55, 316.23, 514.78, 850.00, 1610.12, 355.32, 559.02, 888.84, 1336.27, 671.75, 1185.59, 1971.83, 3196.58, 316.23, 447.21, 670.82, 1024.74, 258.26, 389.31, 598.82, 911.41, 253.24, 382.28, 582.10, 897.67, 315.92, 496.45, 743.21, 1101.34, 2.97, 5.57, 9.33, 19.38, 206.02, 312.74, 458.30, 752.91, 249.52, 374.70, 574.74, 880.91, "
                + "2273.08, 1580.80, 1207.77, 1020.56, 904.92, 781.17, 710.75, 652.16, 607.87, 570.01, 538.84, 510.27, 499.65, 452.46, 447.01, "
                + "3077.30, 2016.40, 1631.50, 1352.36, 1187.27, 1063.01, 971.93, 895.06, 832.00, 789.70, "
                + "457.91, 309.01, 240.99, 209.45, 184.91, 162.64, 148.11, 136.77, 129.51, 126.99, 109.25, 107.76, 104.73, 101.60, 91.6";
        tmpsplits = tmps.split(",");
        for (int t = 0; t < 124; t++) {
            targets[t] = Math.round(Double.valueOf(tmpsplits[t]) * 100.0) / 100.0;
            targets_pbs[t] = targets_grasppr[t] = targets[t]; //bacause they are mostly the same- the exceptions will be amended shortly after this loop
        }        
        targets_pbs[116] = targets_grasppr[116] = 136.80;
        targets_pbs[121] = targets_grasppr[121] = 107.75;
        targets_grasppr[123] = 92.44;

        targets_pbs[48] = 355.52;
        targets_pbs[73] = 5.97;
        targets_pbs[90] = 710.76;
        targets_pbs[92] = 607.88;
        targets_pbs[95] = 510.28;
        targets_pbs[112] = 209.46;
        targets_pbs[114] = 162.65;
        targets_pbs[117] = 129.54;
        targets_pbs[118] = 127.01;
        targets_pbs[120] = 107.78;
        targets_pbs[122] = 101.61;
        targets_pbs[123] = 101.60;
    }//fillin_files_and_ns_and_ps_and_targets
    //----------

    static void test_data(boolean[] x) {
        int k, i, v, minpos, minpos2, minfa, minfa2;
        double mind, mind2;

        boolean[] zx = new boolean[n];
        int[] zP = new int[n], zQ = new int[n], zQinverse = new int[n];
        double zf;
        int[] zC = new int[n];
        int znC;

        int[][] zZ = new int[n][n];
        int[] znZ = new int[n];
        int[] zF = new int[n];
        int[] zFinx = new int[n]; //to indicate the index of each client in  Z

        int[][] zZ2 = new int[n][n];
        int[] znZ2 = new int[n];
        int[] zF2 = new int[n];
        int[] zF2inx = new int[n]; //to indicate the index of each client in  Z2

        double[] zm_p = new double[n];
        int[][] zm_p_v = new int[n][n];
        int[] zm_p_v_cnt = new int[n];

        for (i = 0; i < n; i++) {
            zx[i] = x[i];
        }

        int zp, zq;
        for (v = 0, zp = zq = 0; v < n; v++) {
            if (zx[v]) {
                zP[zp] = v;
                zQinverse[v] = -1;
                zp++;
            } else {
                zQ[zq] = v;
                zQinverse[v] = zq;
                zq++;
            }
        }
        //test p, q, P, Q, Qinverse
        int q = n - p;
        if (zp != p || zq != q || zp + zq != n) {
            ee("WRONG p and/or q");
        }

        //P, Q and Qinverse
        if (!compare_two_lists(zP, P, zp)) {
            ee("WRONG P");
        }
        if (!compare_two_lists(zQ, Q, zq)) {
            ee("WRONG Q");
        }
        if (!compare_two_lists(zQinverse, Qinverse, n)) {
            ee("WRONG Qinverse");
        }

        znC = -1;
        zf = -1;
        for (i = 0; i < n; i++) { //n not p (i.e., for all potential fa)
            znZ[i] = znZ2[i] = 0;
            zm_p[i] = -1;
            zm_p_v_cnt[i] = 0;
        }
        for (v = 0; v < n; v++) {
            minpos = nextpos(v, -1);
            minfa = N[v][minpos];
            mind = D[minfa][v];
            minpos2 = nextpos(v, minpos);
            minfa2 = N[v][minpos2];
            mind2 = D[minfa2][v];

            zF[v] = minfa;
            zZ[minfa][znZ[minfa]] = v;
            zFinx[v] = znZ[minfa];//ma shaa Allah
            znZ[minfa]++;

            zF2[v] = minfa2;
            zZ2[minfa2][znZ2[minfa2]] = v;
            zF2inx[v] = znZ2[minfa2];//ma shaa Allah
            znZ2[minfa2]++;

            if (mind > zf) {
                zf = mind;
                zC[0] = v;
                znC = 1;
            } else {
                if (mind == zf) {
                    zC[znC] = v;
                    znC++;
                }
            }

            if (zx[minfa]) {
                if (mind2 > zm_p[minfa]) {
                    zm_p[minfa] = mind2;
                    zm_p_v[minfa][0] = v;
                    zm_p_v_cnt[minfa] = 1;
                } else {
                    if (mind2 == zm_p[minfa]) {
                        zm_p_v[minfa][zm_p_v_cnt[minfa]] = v;
                        zm_p_v_cnt[minfa]++;
                    }
                }
            }
        }

        if (zf != f) {
            ee("WRONG f");
        }
        if (znC != nC) {
            ee("WRONG nC");
        }
        if (!compare_two_lists(zC, C, znC)) {
            ee("WRONG C");
        }

        for (v = 0; v < n; v++) //F
        {
            if (zF[v] != F[v]) {
                ee("WRONG F[");
            }

            //Finx and Z
            if (znZ[v] != nZ[v]) {
                ee("WRONG nZ");
            }
            if (!compare_two_lists(zZ[v], Z[v], znZ[v])) {
                ee("WRONG Z");
            }
            if (/*zFinx[v] != Finx[v] ||*/zZ[F[v]][zFinx[v]] != v || Z[F[v]][Finx[v]] != v) {
                ee("WRONG Finx or Z");
            }

            //F2
            if (zF2[v] != F2[v]) {
                ee("WRONG F2[");
            }
            //F2inx and Z2
            if (znZ2[v] != nZ2[v]) {
                ee("WRONG nZ2");
            }
            if (!compare_two_lists(zZ2[v], Z2[v], znZ2[v])) {
                ee("WRONG Z2");
            }
            if (/*zF2inx[v] != F2inx[v] ||*/zZ2[zF2[v]][zF2inx[v]] != v || Z2[F2[v]][F2inx[v]] != v) {
                ee("WRONG F2inx or Z2");
            }

            //m_p stuff
            if (zm_p[v] != m_p[v]) {
                ee("WRONG m_p");
            }
            if (zm_p_v_cnt[v] != m_p_v_cnt[v]) {
                ee("WRONG m_p_v_cnt");
            }
            if (!compare_two_lists(zm_p_v[v], m_p_v[v], zm_p_v_cnt[v])) {
                ee("WRONG m_p_v");
            }

            //extra, may not be needed really:
            double tmpmaxF2 = -1;
            for (int k1 = 0; k1 < nZ[v]; k1++) {
                int v1 = Z[v][k1];
                if (D[v1][F2[v1]] > tmpmaxF2) {
                    tmpmaxF2 = D[v1][F2[v1]];
                }
            }
            if (tmpmaxF2 != m_p[v]) {
                ee("inconsistency 8");
            }

            for (int k2 = 0; k2 < m_p_v_cnt[v]; k2++) {
                int v2 = m_p_v[v][k2];
                if (D[v2][F2[v2]] != m_p[v]) {
                    ee("inconsistency 9");
                }
            }

        }//for v
    }//test_data
    //-------

    static boolean compare_two_lists(int[] list1, int[] list2, int cnt) {
        int[] tmplist1 = new int[cnt], tmplist2 = new int[cnt];
        for (int i = 0; i < cnt; i++) {
            tmplist1[i] = list1[i];
            tmplist2[i] = list2[i];
        }
        Arrays.sort(tmplist1);
        Arrays.sort(tmplist2);
        for (int k = 0; k < cnt; k++) {
            if (tmplist1[k] != tmplist2[k]) {
                return false;
            }
        }
        return true;
    }//compare_two_lists (for int[])
    //-------

    static boolean compare_two_lists(double[] list1, double[] list2, int cnt) {
        double[] tmplist1 = new double[cnt], tmplist2 = new double[cnt];
        for (int i = 0; i < cnt; i++) {
            tmplist1[i] = list1[i];
            tmplist2[i] = list2[i];
        }
        Arrays.sort(tmplist1);
        Arrays.sort(tmplist2);
        for (int k = 0; k < cnt; k++) {
            if (tmplist1[k] != tmplist2[k]) {
                return false;
            }
        }
        return true;
    }//compare_two_lists (for double[])
    //----------

}//end class
