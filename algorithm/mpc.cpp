#include "mpc.h"
#include "useful_math.h"

MPC::MPC(double dtIn):QP(nu*mpc_U_N, nc*mpc_U_N) {                          //36，60
    m = 10.0;                                                  //机器人重量
    g = -9.8;                                                  //重力加速度
    miu = 0.5;                                                 //摩擦力系数                           

    max[0] = 1000.0;  max[1] = 1000.0; max[2] = -3.0*m*g;
    max[3] = 1000.0;  max[4] = 1000.0; max[5] = -3.0*m*g;

    min[0] = -1000.0;  min[1] = -1000.0; min[2] = 0.0;
    min[3] = -1000.0;  min[4] = -1000.0; min[5] = 0.0;

    //single rigid body model
    for (int i = 0; i < (mpc_N); i ++){                         //初始化状态变量系数矩阵与控制输入量系数矩阵
        Ac[i].setZero();                                       
        Bc[i].setZero();
        A[i].setZero();
        B[i].setZero();
    }

    Aqp.setZero();
    Aqp1.setZero();
    Bqp1.setZero();
    Bqp.setZero();

    Ufe.setZero();
    Ufe_pre.setZero();

    Xd.setZero();
    X_cur.setZero();
    X_cal.setZero();
    dX_cal.setZero();

    Q = Eigen::MatrixXd::Zero(nx*mpc_N, nx*mpc_N);                             //130， 130
    R.setZero(); M.setZero();
    alpha = 0.0;
    H.setZero();
    c.setZero();

    u_low.setZero(); u_up.setZero();
    As.setZero();
    bs.setZero();

    pCoM.setZero();
    pf2com.setZero(); pe.setZero();
    R_cur.setZero();
    R_w2f.setZero(); R_f2w.setZero();

    qpOASES::Options  option;
    option.printLevel = qpOASES::PL_LOW;
    QP.setOptions(option);

    dt = dtIn;
}

void MPC::set_weight(double u_weight, Eigen::MatrixXd Q_diag, Eigen::MatrixXd R_diag) {
    Eigen::MatrixXd   Q_diag_N = Eigen::MatrixXd::Zero(1, nx*mpc_N);                        //65
    Eigen::MatrixXd   R_diag_N = Eigen::MatrixXd::Zero(1, nu*mpc_U_N);                      //36

    Q = Eigen::MatrixXd::Zero(nx*mpc_N, nx*mpc_N);                  //误差权重   130，130
    R = Eigen::MatrixXd::Zero(nu*mpc_U_N, nu*mpc_U_N);              //系统输入权重 36，36

    alpha = u_weight;                                               //系统输入权重
    
    // 完善K、L矩阵（Q，R）
    for (int i = 0; i < mpc_N; i++) {
        Q_diag_N.block<1,nx>(0, i*nx) = Q_diag;
    }

    for (int i = 0; i < mpc_U_N; i++) {
        R_diag_N.block<1,nu>(0, i*nu) = R_diag;
    }

    for (int i = 0; i < nx*mpc_N; i++) {
        Q(i,i) = Q_diag_N(0,i);
    }

    for (int i = 0; i < nu*mpc_U_N; i++) {
        R(i,i) = R_diag_N(0,i);
    }

	// for (int i = 0; i < mpc_N; i++){
	// 	Q.block<3,3>(i*nx + 3,i*nx + 3) = R_curz[i]*Q.block<3,3>(i*nx + 3,i*nx + 3)*R_curz[i].transpose();
	// 	Q.block<3,3>(i*nx + 6,i*nx + 6) = R_curz[i]*Q.block<3,3>(i*nx + 6,i*nx + 6)*R_curz[i].transpose();
	// 	Q.block<3,3>(i*nx + 9,i*nx + 9) = R_curz[i]*Q.block<3,3>(i*nx + 9,i*nx + 9)*R_curz[i].transpose();
	// }

    // for (int i = 0; i < mpc_U_N; i++){
    //     R.block<3,3>(i*nu,i*nu) = R_curz[i]*R.block<3,3>(i*nu,i*nu)*R_curz[i].transpose();
    //     R.block<3,3>(i*nu + 3,i*nu + 3) = R_curz[i]*R.block<3,3>(i*nu + 3,i*nu + 3)*R_curz[i].transpose();
    //     R.block<3,3>(i*nu + 6,i*nu + 6) = R_curz[i]*R.block<3,3>(i*nu + 6,i*nu + 6)*R_curz[i].transpose();
    //     R.block<3,3>(i*nu + 9,i*nu + 9) = R_curz[i]*R.block<3,3>(i*nu + 9,i*nu + 9)*R_curz[i].transpose();
    // }

}

void MPC::dataBusRead(DataBus &Data) {
    //读取当前状态变量
    X_cur.block<3,1>(0,0) = Data.base_rpy;                  //旋转角
    X_cur.block<3,1>(3,0) = Data.q.block<3,1>(0,0);         //位置
    X_cur.block<3,1>(6,0) = Data.dq.block<3,1>(3,0);        //角速度
    X_cur.block<3,1>(9,0) = Data.dq.block<3,1>(0,0);        //速度   
       
    if (EN) {
        //set Xd 设置期望期望值
        for (int i = 0; i < (mpc_N - 1); i++)
            Xd.block<nx, 1>(nx * i, 0) = Xd.block<nx, 1>(nx * (i + 1), 0);
        for (int j = 0; j < 3; j++)
            Xd(nx * (mpc_N - 1) + j) = Data.js_eul_des(j);
        for (int j = 0; j < 3; j++)
            Xd(nx * (mpc_N - 1) + 3 + j) = Data.js_pos_des(j);
        for (int j = 0; j < 3; j++)
            Xd(nx * (mpc_N - 1) + 6 + j) = Data.js_omega_des(j);
        for (int j = 0; j < 3; j++)
            Xd(nx * (mpc_N - 1) + 9 + j) = Data.js_vel_des(j);
    }
    else{
        for (int j = 0; j < 3; j++)
            Data.js_eul_des(j) = Data.base_rpy(j);
        for (int j = 0; j < 3; j++)
            Data.js_pos_des(j) = Data.base_pos(j);
        for (int j = 0; j < 3; j++)
            Data.js_omega_des(j) = Data.base_omega_W(j);
        for (int j = 0; j < 3; j++)
            Data.js_vel_des(j) = Data.dq(j);

        for (int i = 0; i < mpc_N; i++){
            for (int j = 0; j < 3; j++)
                Xd(nx * i + j) = Data.js_eul_des(j);
            for (int j = 0; j < 3; j++)
                Xd(nx * i + 3 + j) = Data.js_pos_des(j);
            for (int j = 0; j < 3; j++)
                Xd(nx * i + 6 + j) = Data.js_omega_des(j);
            for (int j = 0; j < 3; j++)
                Xd(nx * i + 9 + j) = Data.js_vel_des(j);
        }
    }

    R_cur = eul2Rot(X_cur(0), X_cur(1), X_cur(2));//Data.base_rot;
    for (int i = 0; i < mpc_N; i++) {
        R_curz[i] = Rz3(X_cur(2));
    }

    pCoM = X_cur.block<3,1>(3,0);                   //当前机身位置
    pe.block<3,1>(0,0) = Data.fe_fl_pos_W;           //世界坐标系下左前脚位置
    pe.block<3,1>(3,0) = Data.fe_fr_pos_W;           //世界坐标系下右前脚位置
    pe.block<3,1>(6,0) = Data.fe_rl_pos_W;           //世界坐标系下右后脚位置
    pe.block<3,1>(9,0) = Data.fe_rr_pos_W;           //世界坐标系下左后脚位置
    
    pf2com.block<3,1>(0,0) = pe.block<3,1>(0,0) - pCoM;                 //机身到左前脚的向量表示
    pf2com.block<3,1>(3,0) = pe.block<3,1>(3,0) - pCoM;                 //机身到右前脚的向量表示
    pf2com.block<3,1>(6,0) = pe.block<3,1>(6,0) - pCoM;                 //机身到右后脚的向量表示
    pf2com.block<3,1>(9,0) = pe.block<3,1>(9,0) - pCoM;                 //机身到左后脚的向量表示

    pf2comd.block<3,1>(0,0) = pe.block<3,1>(0,0) - Xd.block<3,1>(0,0);  //机身到左前脚期望位置的向量表示
    pf2comd.block<3,1>(3,0) = pe.block<3,1>(3,0) - Xd.block<3,1>(3,0);  //机身到右前脚期望位置的向量表示
    pf2comd.block<3,1>(6,0) = pe.block<3,1>(6,0) - Xd.block<3,1>(6,0);  //机身到右后脚期望位置的向量表示
    pf2comd.block<3,1>(9,0) = pe.block<3,1>(9,0) - Xd.block<3,1>(9,0);  //机身到左后脚期望位置的向量表示

    Ic = Data.inertia;
    // Ic <<   12.61,  0, 0.37
    //         ,0,  11.15, 0.01
    //         ,0.37,0.01, 2.15;

    legStateCur = Data.legState;
    legStateNext = Data.legStateNext;
    for (int i = 0; i < mpc_N; i++){
        double aa;
        aa = i*dt/0.4;
        double phip;
        phip = Data.phi + aa;
        if (phip > 1)
            legState[i] = legStateNext;
        else
            legState[i] = legStateCur;
    }

    //仍需修改  算法中没有用到这些参数
    Eigen::Matrix<double, 3, 3>     R_slop;
    R_slop = eul2Rot(Data.slop(0), Data.slop(1), Data.slop(2));
    if (legStateCur == DataBus::FrSt)
        R_f2w = Data.fe_r_rot_W;
    else if (legStateCur == DataBus::FlSt)
        R_f2w = Data.fe_l_rot_W;
    else
        R_f2w = R_slop;
    R_w2f = R_f2w.transpose();
}

void MPC::cal() {
    if (EN) {
        //qp pre
		for (int i = 0; i < mpc_N; i++) {
			Ac[i].block<3, 3>(0, 6) = R_curz[i].transpose();
			Ac[i].block<3, 3>(3, 9) = Eigen::MatrixXd::Identity(3,3);
			A[i] = Eigen::MatrixXd::Identity(nx,nx) + dt * Ac[i];
		}
		for (int i = 0; i < mpc_N; i++) {
			pf2comi[i] = pf2com;
			Eigen::Matrix3d Ic_W_inv;                   //惯性张量从本体系到定向本体系的映射的逆
			Ic_W_inv = (R_curz[i] * Ic * R_curz[i].transpose()).inverse();    
			Bc[i].block<3, 3>(6, 0) = Ic_W_inv * CrossProduct_A(pf2comi[i].block<3, 1>(0, 0));
			Bc[i].block<3, 3>(6, 3) = Ic_W_inv * CrossProduct_A(pf2comi[i].block<3, 1>(3, 0));
			Bc[i].block<3, 3>(6, 6) = Ic_W_inv * CrossProduct_A(pf2comi[i].block<3, 1>(6, 0));
			Bc[i].block<3, 3>(6, 9) = Ic_W_inv * CrossProduct_A(pf2comi[i].block<3, 1>(9, 0));
			Bc[i].block<3, 3>(9, 0) = Eigen::MatrixXd::Identity(3,3)/ m;
            Bc[i].block<3, 3>(9, 3) = Eigen::MatrixXd::Identity(3,3)/ m;
            Bc[i].block<3, 3>(9, 6) = Eigen::MatrixXd::Identity(3,3)/ m;
			Bc[i].block<3, 3>(9, 9) = Eigen::MatrixXd::Identity(3,3)/ m;
            Bc[i]((nx - 1), (nu - 1)) = 1.0 / m;
			B[i] = dt * Bc[i];
		}
        //计算预测阶段递推后的状态变量矩阵
		for (int i = 0; i < mpc_N; i++)
			Aqp.block<nx, nx>(i * nx, 0) = Eigen::MatrixXd::Identity(nx,nx);
		for (int i = 0; i < mpc_N; i++)
			for (int j = 0; j < i + 1; j++)
				Aqp.block<nx, nx>(i * nx, 0) = A[j] * Aqp.block<nx, nx>(i * nx, 0);
        
        //计算预测阶段递推后的控制输入量矩阵
		for (int i = 0; i < mpc_N; i++)
			for (int j = 0; j < i + 1; j++)
				Aqp1.block<nx, nx>(i * nx, j * nx) = Eigen::MatrixXd::Identity(nx,nx);
        for (int i = 1; i < mpc_N; i++)
            for (int j = 0; j < i; j++)
                for (int k = j + 1; k < (i + 1); k++)
                    Aqp1.block<nx, nx>(i * nx, j * nx) = A[k] * Aqp1.block<nx, nx>(i * nx, j * nx);
        /*
          |--
          |  I_13
          |  A[1]            I_13
          |  A[1]*A[2]       A[2]       I_13 
          |  A[1]*A[2]*A[3]  A[2]*A[3]  A[3]  I_13
          |  ...........................................
          |  A[1]*A[2]...A[n-1] .............................
          |__

        */
        for (int i = 0; i < mpc_N; i++)
            Bqp1.block<nx, nu>(i * nx, i * nu) = B[i];
        Eigen::MatrixXd Bqp11 = Eigen::MatrixXd::Zero(nu * mpc_N, nu * mpc_U_N);
        Bqp11.setZero();
        Bqp11.block<nu * mpc_U_N, nu * mpc_U_N>(0, 0) = Eigen::MatrixXd::Identity(nu * mpc_U_N, nu * mpc_U_N);
        for (int i = 0; i < (mpc_N - mpc_U_N); i++)
            Bqp11.block<nu, nu>(nu * mpc_U_N + i * nu, nu * (mpc_U_N - 1)) = Eigen::MatrixXd::Identity(nu, nu);

        Eigen::MatrixXd B_tmp = Eigen::MatrixXd::Zero(nx * mpc_N, nu * mpc_U_N);
        B_tmp = Bqp1 * Bqp11;
        Bqp = Aqp1 * B_tmp;                                         //预测阶段递推后的控制输入量矩阵

        H = 2 * (Bqp.transpose() * Q * Bqp + alpha * R);
        // H = Eigen::MatrixXd::Identity(26,26);
        
        c.transpose() = 2 * Bqp.transpose() * Q * (Aqp * X_cur - Xd);
        // c.setZero();


        /*
          代价函数： J(U) = (X - XD)^T * Q * (X - XD) + U^T * R * U
          展开计算整理后： J(U) = 1/2 * U^T * H * U + U^T * c
          其中 H = 2 * Bqp^T * Q * Bqp + R
              c = 2 * Bqp^T * Q * (Aqp * X - XD)
        */

        //friction constraint   摩擦约束
        Eigen::Matrix<double, ncfr, 3> Asfr11;
        Eigen::Matrix<double, ncf*ncfr, nu> Asfr1;
        Eigen::Matrix<double, nc * mpc_U_N, nu * mpc_U_N> Asfr;
        Asfr1.setZero();
        Asfr.setZero();
        // Asfr11 <<
        //         -1.0,  0.0, miu,
        //          0.0, -1.0, miu,
        //          1.0,  0.0, miu,
        //          0.0,  1.0, miu,
        //          0.0,  0.0, 1.0;

        Asfr11 <<
                -1.0,  0.0, -1.0 / sqrt(2.0) * miu,
                 0.0, -1.0, -1.0 / sqrt(2.0) * miu,
                 1.0,  0.0, -1.0 / sqrt(2.0) * miu,
                 0.0,  1.0, -1.0 / sqrt(2.0) * miu;

        Asfr11 = Asfr11 * R_w2f;

        for (int i = 0; i < ncf; i++)
            Asfr1.block<ncfr, 3>(ncfr * i, i * 3) = Asfr11;
        
        for (int i = 0; i < mpc_U_N; i++)
            Asfr.block<nc, nu>(nc * i, i * nu) = Asfr1;
        
        As.block<nc * mpc_U_N, nu * mpc_U_N>(0, 0) = Asfr;

        // As.setZero();

        //qp
        Eigen::Matrix<double, nu * mpc_U_N, 1> Guess_value;
        Guess_value.setZero();
        for (int i = 0; i < mpc_U_N; i++){
            if (legState[i] == DataBus::DSt) {
                Guess_value(i * nu + 2) = -0.5 * m * g;
                Guess_value(i * nu + 5) = -0.5 * m * g;
                Guess_value(i * nu + 8) = -0.5 * m * g;
                Guess_value(i * nu + 11) = -0.5 * m * g;
                Guess_value(i * nu + 12) = m * g;
                for (int j = 0; j < 6; j++) {
                    u_low(i * nu + j) = min[j];
                    u_low(i * nu + j + 6) = min[j];
                    u_up(i * nu + j) = max[j];
                    u_up(i * nu + j + 6) = max[j];
                }
                u_low(i * nu + 12) = m * g;
                u_up(i * nu + 12) = m * g;
            } else if (legState[i] == DataBus::FlSt) {
                Guess_value(i * nu + 2) = -0.5 * m * g;
                Guess_value(i * nu + 5) = 0.0;
                Guess_value(i * nu + 8) = -0.5 * m * g;
				Guess_value(i * nu + 11) = 0.0;
                Guess_value(i * nu + 12) = m * g;
                for (int j = 0; j < 3; j++) {
                    u_low(i * nu + j) = min[j];
                    u_low(i * nu + j + 3) = 0.0;
                    u_low(i * nu + j + 6) = min[j];
                    u_low(i * nu + j + 9) = 0.0;
                    u_up(i * nu + j) = max[j];
                    u_up(i * nu + j + 3) = 0.0;
                    u_up(i * nu + j + 6) = max[j];
                    u_up(i * nu + j + 9) = 0.0;
                }
                u_low(i * nu + 12) = m * g;
                u_up(i * nu + 12) = m * g;
            } else if (legState[i] == DataBus::FrSt) {
				Guess_value(i * nu + 2) = 0.0;
                Guess_value(i * nu + 5) = -0.5 * m * g;
                Guess_value(i * nu + 8) = 0.0;
                Guess_value(i * nu + 11) = -0.5 * m * g;
                Guess_value(i * nu + 12) = m * g;
                for (int j = 0; j < 3; j++) {
                    u_low(i * nu + j) = 0.0;
                    u_low(i * nu + j + 3) = min[j];
                    u_low(i * nu + j + 6) = 0.0;
                    u_low(i * nu + j + 9) = min[j];
                    u_up(i * nu + j) = 0.0;
                    u_up(i * nu + j + 3) = max[j];
                    u_up(i * nu + j + 6) = 0.0;
                    u_up(i * nu + j + 9) = max[j];
                }
                u_low(i * nu + 12) = m * g;
                u_up(i * nu + 12) = m * g;
            }
        }

        qpOASES::returnValue res;
        nWSR = 1000000;
        cpu_time = dt;

        Eigen::Matrix<double, nc * mpc_U_N, 1> lbA, ubA, one_ch_1;
        one_ch_1.setOnes();
        lbA = -1e7 * one_ch_1;
        ubA = 1e7 * one_ch_1;

        // for (int i = 0; i < mpc_U_N; i++){
        //     if(legState[i] == DataBus::DSt){
        //         ubA((nc * i + 4),0) = max[2];
        //         ubA((nc * i + 9),0) = max[2];
        //         ubA((nc * i + 14),0) = max[2];
        //         ubA((nc * i + 19),0) = max[2];
        //     }
        //     else if (legState[i] == DataBus::FlSt) {
        //         ubA((nc * i + 4),0) = max[2];
        //         ubA.block<ncfr,1>(nc * i + 5, 0).setZero();
        //         ubA((nc * i + 14),0) = max[2];
        //         ubA.block<ncfr,1>(nc * i + 15, 0).setZero();
        //     }
        //     else if (legState[i] == DataBus::FrSt) {
        //         ubA.block<ncfr,1>(nc * i, 0).setZero();
        //         ubA((nc * i + 9),0) = max[2];
        //         ubA.block<ncfr,1>(nc * i + 10, 0).setZero();
        //         ubA((nc * i + 19),0) = max[2];
        //     }
        // }

        for (int i = 0; i < mpc_U_N; i++){
            if(legState[i] == DataBus::DSt){
                // ubA((nc * i + 4),0) = max[2];
                // ubA((nc * i + 9),0) = max[2];
                // ubA((nc * i + 14),0) = max[2];
                // ubA((nc * i + 19),0) = max[2];
            }
            else if (legState[i] == DataBus::FlSt) {
                // ubA((nc * i + 4),0) = max[2];
                ubA.block<ncfr,1>(nc * i + 4, 0).setZero();
                // ubA((nc * i + 14),0) = max[2];
                ubA.block<ncfr,1>(nc * i + 12, 0).setZero();
            }
            else if (legState[i] == DataBus::FrSt) {
                ubA.block<ncfr,1>(nc * i, 0).setZero();
                // ubA((nc * i + 9),0) = max[2];
                ubA.block<ncfr,1>(nc * i + 8, 0).setZero();
                // ubA((nc * i + 19),0) = max[2];
            }
        }   

        copy_Eigen_to_real_t(qp_H, H, nu * mpc_U_N, nu * mpc_U_N);
        copy_Eigen_to_real_t(qp_c, c.transpose(), nu * mpc_U_N, 1);
        copy_Eigen_to_real_t(qp_As, As, nc * mpc_U_N, nu * mpc_U_N);
        copy_Eigen_to_real_t(qp_lbA, lbA, nc * mpc_U_N, 1);
        copy_Eigen_to_real_t(qp_ubA, ubA, nc * mpc_U_N, 1);
        copy_Eigen_to_real_t(qp_lu, u_low, nu * mpc_U_N, 1);
        copy_Eigen_to_real_t(qp_uu, u_up, nu * mpc_U_N, 1);
        copy_Eigen_to_real_t(xOpt_iniGuess, Guess_value, nu * mpc_U_N, 1);
        res = QP.init(qp_H, qp_c, qp_As, qp_lu, qp_uu, qp_lbA, qp_ubA, nWSR, &cpu_time, xOpt_iniGuess);
        // res = QP.init(qp_H, qp_c, qp_As, NULL, NULL, qp_lbA, qp_ubA, nWSR, &cpu_time, xOpt_iniGuess);


        qp_Status = qpOASES::getSimpleStatus(res);
        qp_nWSR = nWSR;
        qp_cpuTime = cpu_time;

		if (res!=qpOASES::SUCCESSFUL_RETURN)
		{
			printf("failed!!!!!!!!!!!!!\n");
		}

        qpOASES::real_t xOpt[nu * mpc_U_N];
        QP.getPrimalSolution(xOpt);
        if (qp_Status == 0) {
            for (int i = 0; i < nu * mpc_U_N; i++)
                Ufe(i) = xOpt[i];
        }

        dX_cal = Ac[0] * X_cur + Bc[0] * Ufe.block<nu,1>(0,0);
        Eigen::Matrix<double, nx, 1>    delta_X;
        delta_X.setZero();
        for (int i = 0; i < 3; i++){
            delta_X(i) = 0.5*dX_cal(i+6)*dt*dt;
            delta_X(i+3) = 0.5*dX_cal(i+9)*dt*dt;
            delta_X(i+6) = dX_cal(i+6)*dt;
            delta_X(i+9) = dX_cal(i+9)*dt;
        }

        X_cal = (Aqp * X_cur + Bqp * Ufe).block<nx,1>(nx*0,0) + delta_X;

        Ufe_pre = Ufe.block<nu, 1>(0, 0);
        QP.reset();
        std::cout<<Ufe_pre.transpose()<<"force"<<std::endl;
        // while(1){}
    }
}

void MPC::dataBusWrite(DataBus &Data) {
    Data.Xd = Xd;
    Data.X_cur = X_cur;
    Data.fe_react_tau_cmd = Ufe;
    Data.X_cal = X_cal;
    Data.dX_cal = dX_cal;

    Data.qp_nWSR_MPC = nWSR;
    Data.qp_cpuTime_MPC = cpu_time;
    Data.qpStatus_MPC = qp_Status;

    Data.Fr_ff = Ufe.block<12, 1>(0, 0);

    double k = 5;
    Data.des_ddq.block<2, 1>(0, 0) << dX_cal(9), dX_cal(10);

    Data.des_ddq(5) = k * (Xd(6 + 2) - Data.dq(5));

    Data.des_dq.block<3, 1>(0, 0) << Xd(9 + 0), Xd(9 + 1), Xd(9 + 2);
    Data.des_dq.block<2, 1>(3, 0) << 0.0, 0.0;
    Data.des_dq(5) = Xd(6 + 2);

    Data.des_delta_q.block<2, 1>(0, 0) = Data.des_dq.block<2, 1>(0, 0) * dt;
    Data.des_delta_q(5) = Data.des_dq(5) * dt;

    Data.base_rpy_des << 0.005, 0.00, Xd(2);
    Data.base_pos_des << Xd(3 + 0), Xd(3 + 1), Xd(3 + 2);
}

void    MPC::enable(){
    EN = true;
}
void    MPC::disable(){
    EN = false;
}

bool    MPC::get_ENA(){
    return EN;
}

void MPC::copy_Eigen_to_real_t(qpOASES::real_t* target, Eigen::MatrixXd source, int nRows, int nCols) {
    int count = 0;

    // Strange Behavior: Eigen matrix matrix(count) is stored by columns (not rows)
    // real_t is stored by rows, same to C array
    for (int i = 0; i < nRows; i++) {
        for (int j = 0; j < nCols; j++) {
            target[count] = source(i, j);
            count++;
        }
    }
}