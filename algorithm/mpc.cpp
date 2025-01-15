/*
This is part of OpenLoong Dynamics Control, an open project for the control of biped robot,
Copyright (C) 2024 Humanoid Robot (Shanghai) Co., Ltd, under Apache 2.0.
Feel free to use in any purpose, and cite OpenLoong-Dynamics-Control in any style, to contribute to the advancement of the community.
 <https://gitee.com/panda_23/openloong-dyn-control.git>
 <web@openloong.org.cn>
*/
#include "mpc.h"
#include "useful_math.h"


MPC::MPC(double dtIn):QP(nu*ch, nc*ch) {
    m = 28.0;
    g = -9.8;
    miu = 0.5;

    max[0] = -0.1*m*g;  max[1] = -0.1*m*g; max[2] = -1.0*m*g;
    max[3] = -0.1*m*g;  max[4] = -0.1*m*g; max[5] = -1.0*m*g;

    min[0] = 0.1*m*g;  min[1] = 0.1*m*g; min[2] = 0.0;
    min[3] = 0.1*m*g;  min[4] = 0.1*m*g; min[5] = 0.0;


    //single rigid body model
    for (int i = 0; i < (mpc_N); i ++){
        Ac[i].setZero();
        Bc[i].setZero();
        A[i].setZero();
        B[i].setZero();
    }
    Cc.setZero();
    C.setZero();

    Aqp.setZero();
    Aqp1.setZero();
    Bqp1.setZero();
    Bqp.setZero();
    Cqp1.setZero();
    Cqp.setZero();

    Ufe.setZero();
    Ufe_pre.setZero();

    Xd.setZero();
    X_cur.setZero();
    X_cal.setZero();
    dX_cal.setZero();

    L = Eigen::MatrixXd::Zero(nx*mpc_N, nx*mpc_N);
    K.setZero(); M.setZero();
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

void MPC::set_weight(double u_weight, Eigen::MatrixXd L_diag, Eigen::MatrixXd K_diag) {
    Eigen::MatrixXd   L_diag_N = Eigen::MatrixXd::Zero(1, nx*mpc_N);
    Eigen::MatrixXd   K_diag_N = Eigen::MatrixXd::Zero(1, nu*ch);

    L = Eigen::MatrixXd::Zero(nx*mpc_N, nx*mpc_N);
    K = Eigen::MatrixXd::Zero(nu*ch, nu*ch);

    alpha = u_weight;
    for (int i = 0; i < mpc_N; i++) {
        L_diag_N.block<1,nx>(0, i*nx) = L_diag;
    }

    for (int i = 0; i < ch; i++) {
        K_diag_N.block<1,nu>(0, i*nu) = K_diag;
    }

    for (int i = 0; i < nx*mpc_N; i++) {
        L(i,i) = L_diag_N(0,i);
    }

    for (int i = 0; i < nu*ch; i++) {
        K(i,i) = K_diag_N(0,i);
    }

	for (int i = 0; i < mpc_N; i++){
		L.block<3,3>(i*nx + 3,i*nx + 3) = R_curz[i]*L.block<3,3>(i*nx + 3,i*nx + 3)*R_curz[i].transpose();
		L.block<3,3>(i*nx + 6,i*nx + 6) = R_curz[i]*L.block<3,3>(i*nx + 6,i*nx + 6)*R_curz[i].transpose();
		L.block<3,3>(i*nx + 9,i*nx + 9) = R_curz[i]*L.block<3,3>(i*nx + 9,i*nx + 9)*R_curz[i].transpose();
	}

    for (int i = 0; i < ch; i++){
        K.block<3,3>(i*nu,i*nu) = R_curz[i]*K.block<3,3>(i*nu,i*nu)*R_curz[i].transpose();
        K.block<3,3>(i*nu + 3,i*nu + 3) = R_curz[i]*K.block<3,3>(i*nu + 3,i*nu + 3)*R_curz[i].transpose();
        K.block<3,3>(i*nu + 6,i*nu + 6) = R_curz[i]*K.block<3,3>(i*nu + 6,i*nu + 6)*R_curz[i].transpose();
        K.block<3,3>(i*nu + 9,i*nu + 9) = R_curz[i]*K.block<3,3>(i*nu + 9,i*nu + 9)*R_curz[i].transpose();
    }
}


void MPC::dataBusRead(DataBus &Data) {
    //set value
    X_cur.block<3,1>(0,0) = Data.base_rpy;
    X_cur.block<3,1>(3,0) = Data.q.block<3,1>(0,0);
    X_cur.block<3,1>(6,0) = Data.dq.block<3,1>(3,0);
    X_cur.block<3,1>(9,0) = Data.dq.block<3,1>(0,0);
    if (EN) {
        //set Xd
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
        for (int i = 0; i < mpc_N; i++){
            for (int j = 0; j < 3; j++)
                Xd(nx * i + j) = X_cur(j);//Data.js_eul_des(j);
            for (int j = 0; j < 3; j++)
                Xd(nx * i + 3 + j) = X_cur(3 + j);//Data.js_pos_des(j);
            for (int j = 0; j < 3; j++)
                Xd(nx * i + 6 + j) = X_cur(6 + j);//Data.js_omega_des(j);
            for (int j = 0; j < 3; j++)
                Xd(nx * i + 9 + j) = X_cur(9 + j);//Data.js_vel_des(j);
        }
//		for (int j = 0; j < 3; j++)
//			Data.js_eul_des(j) = X_cur(j);//
//		for (int j = 0; j < 3; j++)
//			Data.js_pos_des(j) = X_cur(3 + j);//
//		for (int j = 0; j < 3; j++)
//			Data.js_omega_des(j) = X_cur(6 + j);//;
//		for (int j = 0; j < 3; j++)
//			Data.js_vel_des(j) = X_cur(9 + j);//;
    }
	
    R_cur = eul2Rot(X_cur(0), X_cur(1), X_cur(2));//Data.base_rot;
    for (int i = 0; i < mpc_N; i++) {
        R_curz[i] = Rz3(X_cur(2));
    }
    pCoM = X_cur.block<3,1>(3,0);
    pe.block<3,1>(0,0) = Data.fe_fl_pos_W;       
    pe.block<3,1>(3,0) = Data.fe_fr_pos_W;           
    pe.block<3,1>(6,0) = Data.fe_rl_pos_W;           
    pe.block<3,1>(9,0) = Data.fe_rr_pos_W;     

    pf2com.block<3,1>(0,0) = pe.block<3,1>(0,0) - pCoM;             
    pf2com.block<3,1>(3,0) = pe.block<3,1>(3,0) - pCoM;    
    pf2com.block<3,1>(6,0) = pe.block<3,1>(6,0) - pCoM;  
    pf2com.block<3,1>(9,0) = pe.block<3,1>(9,0) - pCoM;

    pf2comd.block<3,1>(0,0) = pe.block<3,1>(0,0) - Xd.block<3,1>(0,0);  
    pf2comd.block<3,1>(3,0) = pe.block<3,1>(3,0) - Xd.block<3,1>(3,0); 
    pf2comd.block<3,1>(6,0) = pe.block<3,1>(6,0) - Xd.block<3,1>(6,0);  
    pf2comd.block<3,1>(9,0) = pe.block<3,1>(9,0) - Xd.block<3,1>(9,0); 

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

    // Eigen::Matrix<double, 3, 3>     R_slop;
    // R_slop = eul2Rot(Data.slop(0), Data.slop(1), Data.slop(2));
    // if (legStateCur == DataBus::RSt)
    //     R_f2w = Data.fe_r_rot_W;
    // else if (legStateCur == DataBus::LSt)
    //     R_f2w = Data.fe_l_rot_W;
    // else
    //     R_f2w = R_slop;
    // R_w2f = R_f2w.transpose();
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
			Eigen::Matrix3d Ic_W_inv;
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
		for (int i = 0; i < mpc_N; i++)
			Aqp.block<nx, nx>(i * nx, 0) = Eigen::MatrixXd::Identity(nx,nx);
		for (int i = 0; i < mpc_N; i++)
			for (int j = 0; j < i + 1; j++)
				Aqp.block<nx, nx>(i * nx, 0) = A[j] * Aqp.block<nx, nx>(i * nx, 0);

		for (int i = 0; i < mpc_N; i++)
			for (int j = 0; j < i + 1; j++)
				Aqp1.block<nx, nx>(i * nx, j * nx) = Eigen::MatrixXd::Identity(nx,nx);
        for (int i = 1; i < mpc_N; i++)
            for (int j = 0; j < i; j++)
                for (int k = j + 1; k < (i + 1); k++)
                    Aqp1.block<nx, nx>(i * nx, j * nx) = A[k] * Aqp1.block<nx, nx>(i * nx, j * nx);

        for (int i = 0; i < mpc_N; i++)
            Bqp1.block<nx, nu>(i * nx, i * nu) = B[i];
        Eigen::MatrixXd Bqp11 = Eigen::MatrixXd::Zero(nu * mpc_N, nu * ch);
        Bqp11.setZero();
        Bqp11.block<nu * ch, nu * ch>(0, 0) = Eigen::MatrixXd::Identity(nu * ch, nu * ch);
        for (int i = 0; i < (mpc_N - ch); i++)
            Bqp11.block<nu, nu>(nu * ch + i * nu, nu * (ch - 1)) = Eigen::MatrixXd::Identity(nu, nu);

        Eigen::MatrixXd B_tmp = Eigen::MatrixXd::Zero(nx * mpc_N, nu * ch);
        B_tmp = Bqp1 * Bqp11;
        Bqp = Aqp1 * B_tmp;

		Eigen::Matrix<double, nu*ch, 1>		delta_U;
		delta_U.setZero();
		for (int i = 0; i < ch; i++){
			if (legState[i] == DataBus::FlSt){
                delta_U(nu*i + 2) = 0.5*m*g;
                delta_U(nu*i + 11) = 0.5*m*g;
            }
			else if (legState[i] == DataBus::FrSt){
                delta_U(nu*i + 5) = 0.5*m*g;
                delta_U(nu*i + 8) = 0.5*m*g;
            }
				
			else{
				delta_U(nu*i + 2) = 0.25*m*g;
                delta_U(nu*i + 5) = 0.25*m*g;
				delta_U(nu*i + 8) = 0.25*m*g;
                delta_U(nu*i + 11) = 0.25*m*g;
			}
		}

		H = 2 * (Bqp.transpose() * L * Bqp + alpha * K) + 1e-10*Eigen::MatrixXd::Identity(nx*mpc_N, nx*mpc_N);
		c = 2 * Bqp.transpose() * L * (Aqp * X_cur - Xd) + 2 * alpha * K * delta_U;

        //friction constraint
        Eigen::Matrix<double, ncfr, 3> Asfr11;
        Eigen::Matrix<double, ncf*ncfr, nu> Asfr1;
        Eigen::Matrix<double, nc * ch, nu * ch> Asfr;
        Asfr1.setZero();
        Asfr.setZero();

        Asfr11 <<
                -1.0,  0.0, miu,
                 0.0, -1.0, miu,
                 1.0,  0.0, miu,
                 0.0,  1.0, miu,
                 0.0,  0.0, 1.0;

        // Asfr11 = Asfr11 * R_w2f;

        for (int i = 0; i < ncf; i++)
            Asfr1.block<ncfr, 3>(ncfr * i, i * 3) = Asfr11;
        
        for (int i = 0; i < ch; i++)
            Asfr.block<nc, nu>(nc * i, i * nu) = Asfr1;
        
        As.block<nc * ch, nu * ch>(0, 0) = Asfr;

        //qp
        Eigen::Matrix<double, nu * ch, 1> Guess_value;
        Guess_value.setZero();
        for (int i = 0; i < ch; i++){
            if (legState[i] == DataBus::DSt) {
                Guess_value(i * nu + 2) = -0.25 * m * g;
                Guess_value(i * nu + 5) = -0.25 * m * g;
                Guess_value(i * nu + 8) = -0.25 * m * g;
                Guess_value(i * nu + 11) = -0.25 * m * g;
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
                // Guess_value(i * nu + 5) = 0.0;
                // Guess_value(i * nu + 8) = 0.0;
                Guess_value(i * nu + 11) = -0.5 * m * g;
                Guess_value(i * nu + 12) = m * g;
                for (int j = 0; j < 3; j++) {
                    u_low(i * nu + j) = min[j];
                    // u_low(i * nu + j + 3) = 0.0;
                    // u_low(i * nu + j + 6) = 0.0;
                    u_low(i * nu + j + 9) = min[j];
                    
                    u_up(i * nu + j) = max[j];
                    // u_up(i * nu + j + 3) = 0.0;
                    // u_up(i * nu + j + 6) = 0.0;
                    u_up(i * nu + j + 9) = max[j];
                }
                u_low(i * nu + 12) = m * g;
                u_up(i * nu + 12) = m * g;
            } else if (legState[i] == DataBus::FrSt) {
				// Guess_value(i * nu + 2) = 0.0;
                Guess_value(i * nu + 5) = -0.5 * m * g;
                Guess_value(i * nu + 8) = -0.5 * m * g;
                // Guess_value(i * nu + 11) = 0.0;
                Guess_value(i * nu + 12) = m * g;
                for (int j = 0; j < 3; j++) {
                    // u_low(i * nu + j) = 0.0;
                    u_low(i * nu + j + 3) = min[j];
                    u_low(i * nu + j + 6) = min[j];
                    // u_low(i * nu + j + 9) = 0.0;
                    
                    // u_up(i * nu + j) = 0.0;
                    u_up(i * nu + j + 3) = max[j];
                    u_up(i * nu + j + 6) = max[j];
                    // u_up(i * nu + j + 9) = 0.0;
                }
                u_low(i * nu + 12) = m * g;
                u_up(i * nu + 12) = m * g;
            }
        }

        qpOASES::returnValue res;
        nWSR = 1000000;
        cpu_time = dt;

        Eigen::Matrix<double, nc * ch, 1> lbA, ubA, one_ch_1;
        one_ch_1.setOnes();
        lbA = -1e7 * one_ch_1;
        ubA = 1e7 * one_ch_1;

        for (int i = 0; i < ch; i++){
            if(legState[i] == DataBus::DSt){
                ubA((nc * i + 4),0) = max[2];
                ubA((nc * i + 9),0) = max[2];
                ubA((nc * i + 14),0) = max[2];
                ubA((nc * i + 19),0) = max[2];
            }
            else if (legState[i] == DataBus::FlSt) {
                ubA((nc * i + 4),0) = max[2];
                ubA.block<ncfr,1>(nc * i + 5, 0).setZero();
                ubA((nc * i + 14),0) = max[2];
                ubA.block<ncfr,1>(nc * i + 15, 0).setZero();
            }
            else if (legState[i] == DataBus::FrSt) {
                ubA.block<ncfr,1>(nc * i, 0).setZero();
                ubA((nc * i + 9),0) = max[2];
                ubA.block<ncfr,1>(nc * i + 10, 0).setZero();
                ubA((nc * i + 19),0) = max[2];
            }
        }

        copy_Eigen_to_real_t(qp_H, H, nu * ch, nu * ch);
        copy_Eigen_to_real_t(qp_c, c, nu * ch, 1);
        copy_Eigen_to_real_t(qp_As, As, nc * ch, nu * ch);
        copy_Eigen_to_real_t(qp_lbA, lbA, nc * ch, 1);
        copy_Eigen_to_real_t(qp_ubA, ubA, nc * ch, 1);
        copy_Eigen_to_real_t(qp_lu, u_low, nu * ch, 1);
        copy_Eigen_to_real_t(qp_uu, u_up, nu * ch, 1);
        copy_Eigen_to_real_t(xOpt_iniGuess, Guess_value, nu * ch, 1);
        res = QP.init(qp_H, qp_c, qp_As, qp_lu, qp_uu, NULL, NULL, nWSR, &cpu_time, xOpt_iniGuess);

        qp_Status = qpOASES::getSimpleStatus(res);
        qp_nWSR = nWSR;
        qp_cpuTime = cpu_time;

		if (res!=qpOASES::SUCCESSFUL_RETURN)
		{
//			printf("failed!!!!!!!!!!!!!\n");
		}

        qpOASES::real_t xOpt[nu * ch];
        QP.getPrimalSolution(xOpt);
        if (qp_Status == 0) {
            for (int i = 0; i < nu * ch; i++)
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

        std::cout<<"F: "<<Ufe_pre<<std::endl;
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

