/*
This is part of OpenLoong Dynamics Control, an open project for the control of biped robot,
Copyright (C) 2024 Humanoid Robot (Shanghai) Co., Ltd, under Apache 2.0.
Feel free to use in any purpose, and cite OpenLoong-Dynamics-Control in any style, to contribute to the advancement of the community.
 <https://atomgit.com/openloong/openloong-dyn-control.git>
 <web@openloong.org.cn>
*/
//
// Created by boxing on 23-12-29.
//

#include "wbc_priority.h"
#include "iostream"

// QP_nvIn=18, QP_ncIn=22
WBC_priority::WBC_priority(int model_nv_In, int QP_nvIn, int QP_ncIn, double miu_In, double dt) : QP_prob(QP_nvIn,
                                                                                                          QP_ncIn) {
    timeStep = dt;
    model_nv = model_nv_In;
    miu = miu_In;
    QP_nc = QP_ncIn;
    QP_nv = QP_nvIn;
    Sf = Eigen::MatrixXd::Zero(6, model_nv);
    Sf.block<6, 6>(0, 0) = Eigen::MatrixXd::Identity(6, 6);
    St_qpV2 = Eigen::MatrixXd::Zero(model_nv, model_nv - 6); // 6 means the dims of floating base
    St_qpV2.block(6, 0, model_nv - 6, model_nv - 6) = Eigen::MatrixXd::Identity(model_nv - 6, model_nv - 6);

    St_qpV1 = Eigen::MatrixXd::Zero(model_nv, 6); // 6 means the dims of delta_b
    St_qpV1.block<6, 6>(0, 0) = Eigen::MatrixXd::Identity(6, 6);

    // defined in body frame
    f_z_low = 10;
    f_z_upp = 300;

    tau_upp_stand_L << 10, 20, 40;    // foot end contact torque limit for stand state, in body frame
    tau_low_stand_L << -10, -20, -40;

    tau_upp_walk_L << 25, 40, 40; // foot end contact torque limit for walk state, in body frame
    tau_low_walk_L << -25, -40, -40;

    qpOASES::Options options;
    options.setToMPC();
    //options.setToReliable();
    options.printLevel = qpOASES::PL_LOW;
    QP_prob.setOptions(options);

    eigen_xOpt = Eigen::VectorXd::Zero(QP_nv);
    eigen_ddq_Opt = Eigen::VectorXd::Zero(model_nv);
    eigen_fr_Opt = Eigen::VectorXd::Zero(12);
    eigen_tau_Opt = Eigen::VectorXd::Zero(model_nv - 6);

    delta_q_final_kin = Eigen::VectorXd::Zero(model_nv);
    dq_final_kin = Eigen::VectorXd::Zero(model_nv);;
    ddq_final_kin = Eigen::VectorXd::Zero(model_nv);

    base_rpy_cur = Eigen::VectorXd::Zero(3);

    //  WBC task defined and order build
    ///------------ walk --------------
    kin_tasks.addTask("static_Contact");
    kin_tasks.addTask("Roll_Pitch_Yaw");
    // kin_tasks_walk.addTask("Roll_Pitch_Yaw_Pz");
    // kin_tasks_walk.addTask("RedundantJoints");
    kin_tasks.addTask("PxPyPz");
    kin_tasks.addTask("SwingLeg");
    // kin_tasks_walk.addTask("HandTrack");
    // kin_tasks_walk.addTask("HandTrackJoints");
    // kin_tasks_walk.addTask("PosRot");

    std::vector<std::string> taskOrder;
    // taskOrder_walk.emplace_back("RedundantJoints");
    taskOrder.emplace_back("static_Contact");
    taskOrder.emplace_back("Roll_Pitch_Yaw");
//    taskOrder_walk.emplace_back("Roll_Pitch_Yaw_Pz");
    taskOrder.emplace_back("PxPyPz");
    // taskOrder_walk.emplace_back("PosRot");
    taskOrder.emplace_back("SwingLeg");
    // taskOrder_walk.emplace_back("HandTrackJoints");


    kin_tasks.buildPriority(taskOrder);

    // ///-------- stand ------------
    // kin_tasks_stand.addTask("static_Contact");
    // kin_tasks_stand.addTask("CoMTrack");
    // kin_tasks_stand.addTask("HandTrackJoints");
    // kin_tasks_stand.addTask("HipRPY");
    // kin_tasks_stand.addTask("HeadRP");
    // kin_tasks_stand.addTask("Pz");
    // kin_tasks_stand.addTask("CoMXY_HipRPY");
    // kin_tasks_stand.addTask("Roll_Pitch_Yaw");
    // kin_tasks_stand.addTask("fixedWaist");


//     std::vector<std::string> taskOrder_stand;

// //    taskOrder_stand.emplace_back("fixedWaist");
//     taskOrder_stand.emplace_back("static_Contact");
// //    taskOrder_stand.emplace_back("CoMTrack");
// //    taskOrder_stand.emplace_back("HipRPY");
//     taskOrder_stand.emplace_back("CoMXY_HipRPY");
//     taskOrder_stand.emplace_back("Pz");
//     taskOrder_stand.emplace_back("HandTrackJoints");
//     taskOrder_stand.emplace_back("HeadRP");

//     kin_tasks_stand.buildPriority(taskOrder_stand);
}

void WBC_priority::dataBusRead(const DataBus &robotState) {
    // foot-end offset posture
    fe_L_rot_L_off = robotState.fe_L_rot_L_off;
    fe_R_rot_L_off = robotState.fe_R_rot_L_off;

    // deisred values
    base_rpy_des = robotState.base_rpy_des;
    base_rpy_cur << robotState.rpy[0], robotState.rpy[1], robotState.rpy[2];
    base_pos_des = robotState.base_pos_des;

    swing_F_fe_pos_des_W = robotState.swing_F_fe_pos_des_W;
    swing_R_fe_pos_des_W = robotState.swing_R_fe_pos_des_W;

    stance_F_fe_pos_cur_W = robotState.stance_F_fe_pos_cur_W;
    stance_R_fe_pos_cur_W = robotState.stance_R_fe_pos_cur_W;

    stance_fe_rot_cur_W = robotState.stance_fe_rot_cur_W;

    des_ddq = robotState.des_ddq;
    des_dq = robotState.des_dq;
    des_delta_q = robotState.des_delta_q;
    des_q = robotState.des_q;

    // state update
    J_base = robotState.J_base;
    dJ_base = robotState.dJ_base;
    base_rot = robotState.base_rot;
    base_pos = robotState.base_pos;

    Jfe = Eigen::MatrixXd::Zero(12, model_nv);
    Jfe.block(0, 0, 3, model_nv) = robotState.J_fl.block(0,0,3,model_nv);
    Jfe.block(3, 0, 3, model_nv) = robotState.J_fr.block(0,0,3,model_nv);
    Jfe.block(6, 0, 3, model_nv) = robotState.J_rl.block(0,0,3,model_nv);
    Jfe.block(9, 0, 3, model_nv) = robotState.J_rr.block(0,0,3,model_nv);
    dJfe = Eigen::MatrixXd::Zero(12, model_nv);
    dJfe.block(0, 0, 3, model_nv) = robotState.dJ_fl.block(0,0,3,model_nv);
    dJfe.block(3, 0, 3, model_nv) = robotState.dJ_fr.block(0,0,3,model_nv);
    dJfe.block(6, 0, 3, model_nv) = robotState.dJ_rl.block(0,0,3,model_nv);
    dJfe.block(9, 0, 3, model_nv) = robotState.dJ_rr.block(0,0,3,model_nv);

    Fr_ff = robotState.Fr_ff;
    dyn_M = robotState.dyn_M;
    dyn_M_inv = robotState.dyn_M_inv;
    dyn_Ag = robotState.dyn_Ag;
    dyn_dAg = robotState.dyn_dAg;
    dyn_Non = robotState.dyn_Non;
    dq = robotState.dq;
    q = robotState.q;
    legStateCur = robotState.legState;
    motionStateCur = robotState.motionState;

    Jc = Eigen::MatrixXd::Zero(6,model_nv);
    dJc = Eigen::MatrixXd::Zero(6,model_nv);   
    Jsw = Eigen::MatrixXd::Zero(6,model_nv);
    dJsw = Eigen::MatrixXd::Zero(6,model_nv);

    if (legStateCur == DataBus::FlSt) {
        Jc.block(0,0,3,model_nv) = robotState.J_fl.block(0,0,3,model_nv);
        Jc.block(3,0,3,model_nv) = robotState.J_rr.block(0,0,3,model_nv);
        dJc.block(0,0,3,model_nv) = robotState.dJ_fl.block(0,0,3,model_nv);
        dJc.block(3,0,3,model_nv) = robotState.dJ_rr.block(0,0,3,model_nv);

        Jsw.block(0,0,3,model_nv) = robotState.J_fr.block(0,0,3,model_nv);
        Jsw.block(3,0,3,model_nv) = robotState.J_rl.block(0,0,3,model_nv);
        dJsw.block(0,0,3,model_nv) = robotState.dJ_fr.block(0,0,3,model_nv);
        dJsw.block(3,0,3,model_nv) = robotState.dJ_rl.block(0,0,3,model_nv);

        F_fe_pos_sw_W = robotState.fe_fr_pos_W;
        R_fe_pos_sw_W = robotState.fe_rl_pos_W;
    } else {
        Jc.block(0,0,3,model_nv) = robotState.J_fr.block(0,0,3,model_nv);
        Jc.block(3,0,3,model_nv) = robotState.J_rl.block(0,0,3,model_nv);
        dJc.block(0,0,3,model_nv) = robotState.dJ_fr.block(0,0,3,model_nv);
        dJc.block(3,0,3,model_nv) = robotState.dJ_rl.block(0,0,3,model_nv);

        Jsw.block(0,0,3,model_nv) = robotState.J_fl.block(0,0,3,model_nv);
        Jsw.block(3,0,3,model_nv) = robotState.J_rr.block(0,0,3,model_nv);
        dJsw.block(0,0,3,model_nv) = robotState.dJ_fl.block(0,0,3,model_nv);
        dJsw.block(3,0,3,model_nv) = robotState.dJ_rr.block(0,0,3,model_nv);

        F_fe_pos_sw_W = robotState.fe_fl_pos_W;
        R_fe_pos_sw_W = robotState.fe_rr_pos_W;
    }

}

void WBC_priority::dataBusWrite(DataBus &robotState) {
    robotState.wbc_ddq_final = eigen_ddq_Opt;
    robotState.wbc_tauJointRes = tauJointRes;
    robotState.wbc_FrRes = eigen_fr_Opt;
    robotState.qp_cpuTime = cpu_time;
    robotState.qp_nWSR = nWSR;
    robotState.qp_status = qpStatus;

    robotState.wbc_delta_q_final = delta_q_final_kin;
    robotState.wbc_dq_final = dq_final_kin;
    robotState.wbc_ddq_final = ddq_final_kin;

    robotState.qp_status = qpStatus;
    robotState.qp_nWSR = nWSR;
    robotState.qp_cpuTime = cpu_time;
}

// QP problem contains joint torque, QP_nv=6+12, QP_nc=22;
void WBC_priority::computeTau() {
    // constust the QP problem, refer to the md file for more details
    Eigen::MatrixXd eigen_qp_A1 = Eigen::MatrixXd::Zero(6, QP_nv);// 18 means the sum of dims of delta_r and delta_Fr
    eigen_qp_A1.block<6, 6>(0, 0) = Sf * dyn_M * St_qpV1;
    eigen_qp_A1.block<6, 12>(0, 6) = -Sf * Jfe.transpose();
    Eigen::VectorXd eqRes = Eigen::VectorXd::Zero(6);

    eqRes = -Sf * dyn_M * ddq_final_kin - Sf * dyn_Non + Sf * Jfe.transpose() * Fr_ff;

//     Eigen::Matrix3d Rfe;
//     if (motionStateCur==DataBus::Stand){
//         Rfe = fe_l_rot_cur_W;
//     } else{
//         Rfe = stance_fe_rot_cur_W;
//     }

//     Eigen::Matrix<double,12,12> Mw2b;
//     Mw2b.setZero();
//     Mw2b.block(0,0,3,3)=Rfe.transpose();
//     Mw2b.block(3,3,3,3)=Rfe.transpose();
//     Mw2b.block(6,6,3,3)=Rfe.transpose();
//     Mw2b.block(9,9,3,3)=Rfe.transpose();

//     Eigen::MatrixXd W = Eigen::MatrixXd::Zero(16, 12);
//     W(0, 0) = 1;
//     W(0, 2) = sqrt(2) / 2.0 * miu;
//     W(1, 0) = -1;
//     W(1, 2) = sqrt(2) / 2.0 * miu;
//     W(2, 1) = 1;
//     W(2, 2) = sqrt(2) / 2.0 * miu;
//     W(3, 1) = -1;
//     W(3, 2) = sqrt(2) / 2.0 * miu;
//     W.block<4, 4>(4, 2) = Eigen::MatrixXd::Identity(4, 4);
//     W.block<8, 6>(8, 6) = W.block<8, 6>(0, 0);
//     W=W*Mw2b;

//     Eigen::VectorXd f_low = Eigen::VectorXd::Zero(16);
//     Eigen::VectorXd f_upp = Eigen::VectorXd::Zero(16);
//     Eigen::Vector3d tau_upp_fe, tau_low_fe;
//     if (motionStateCur==DataBus::Stand) {
//         tau_upp_fe = tau_upp_stand_L;
//         tau_low_fe = tau_low_stand_L;
//     }
//     else
//     {
//         tau_upp_fe = tau_upp_walk_L;
//         tau_low_fe = tau_low_walk_L;
//     }
// //    std::cout<<"wbc_computeTau, st_fe_rot"<<std::endl<<stance_fe_rot_cur_W<<std::endl;

//     f_upp.block<8, 1>(0, 0) << 1e10, 1e10, 1e10, 1e10,
//             f_z_upp, tau_upp_fe(0), tau_upp_fe(1), tau_upp_fe(2);
//     f_upp.block<8, 1>(8, 0) = f_upp.block<8, 1>(0, 0);
//     f_low.block<8, 1>(0, 0) << 0, 0, 0, 0,
//             f_z_low, tau_low_fe(0), tau_low_fe(1), tau_low_fe(2);
//     f_low.block<8, 1>(8, 0) = f_low.block<8, 1>(0, 0);

//     if (motionStateCur == DataBus::Walk  || motionStateCur==DataBus::Walk2Stand) {
//         if (legStateCur == DataBus::LSt) {
//             f_upp(12) = 0;
//             f_upp(13) = 0;
//             f_upp(14) = 0;
//             f_upp(15) = 0;

//             f_low(12) = 0;
//             f_low(13) = 0;
//             f_low(14) = 0;
//             f_low(15) = 0;

//             f_low(8) = -1e-7;
//             f_low(9) = -1e-7;
//             f_low(10) = -1e-7;
//             f_low(11) = -1e-7;
//         } else if (legStateCur == DataBus::RSt) {
//             f_upp(4) = 0;
//             f_upp(5) = 0;
//             f_upp(6) = 0;
//             f_upp(7) = 0;

//             f_low(4) = 0;
//             f_low(5) = 0;
//             f_low(6) = 0;
//             f_low(7) = 0;

//             f_low(0) = -1e-7;
//             f_low(1) = -1e-7;
//             f_low(2) = -1e-7;
//             f_low(3) = -1e-7;
//         }
//     }

//     Eigen::MatrixXd eigen_qp_A2 = Eigen::MatrixXd::Zero(16, 18);
//     eigen_qp_A2.block<16, 12>(0, 6) = W;
//     Eigen::VectorXd neqRes_low = Eigen::VectorXd::Zero(16);
//     Eigen::VectorXd neqRes_upp = Eigen::VectorXd::Zero(16);
    Eigen::MatrixXd W=Eigen::MatrixXd::Zero(20,12);
    W(0,0)=1;W(0,2)= sqrt(2)/2.0*miu;
    W(1,0)=-1;W(1,2)= sqrt(2)/2.0*miu;
    W(2,1)=1;W(2,2)= sqrt(2)/2.0*miu;
    W(3,1)=-1;W(3,2)= sqrt(2)/2.0*miu;
    W(4,2)=1;
    W.block<5,3>(5,3)=W.block<5,3>(0,0);
    W.block<10,6>(10,6)=W.block<10,6>(0,0);
    Eigen::VectorXd f_low=Eigen::VectorXd::Zero(20);
    Eigen::VectorXd f_upp=Eigen::VectorXd::Zero(20);

    f_upp.block<5,1>(0,0)<<1e10,1e10,1e10,1e10,f_z_upp;
    f_upp.block<5,1>(5,0)=f_upp.block<5,1>(0,0);
    f_upp.block<10,1>(10,0)=f_upp.block<10,1>(0,0);

    f_low.block<5,1>(0,0)<<0,0,0,0,f_z_low;
    f_low.block<5,1>(5,0)=f_low.block<5,1>(0,0);
    f_low.block<10,1>(10,0)=f_low.block<10,1>(0,0);

    if (legStateCur==DataBus::FlSt)
    {
        f_upp(9)=0;
        f_upp(14)=0;
        f_low(9)=0;
        f_low(14)=0;
    }
    else if (legStateCur==DataBus::FrSt)
    {
        f_upp(4)=0;
        f_upp(19)=0;
        f_low(4)=0;
        f_low(19)=0;
    }

    Eigen::MatrixXd eigen_qp_A2=Eigen::MatrixXd::Zero(20, 18);
    eigen_qp_A2.block<20,12>(0,6)=W;
    Eigen::VectorXd neqRes_low=Eigen::VectorXd::Zero(20);
    Eigen::VectorXd neqRes_upp=Eigen::VectorXd::Zero(20);

    neqRes_low = f_low - W * Fr_ff;
    neqRes_upp = f_upp - W * Fr_ff;

    Eigen::MatrixXd eigen_qp_A_final = Eigen::MatrixXd::Zero(QP_nc, QP_nv);
    eigen_qp_A_final.block<6, 18>(0, 0) = eigen_qp_A1;
    eigen_qp_A_final.block<20, 18>(6, 0) = eigen_qp_A2;

    Eigen::VectorXd eigen_qp_lbA = Eigen::VectorXd::Zero(26);
    Eigen::VectorXd eigen_qp_ubA = Eigen::VectorXd::Zero(26);

    eigen_qp_lbA.block<6, 1>(0, 0) = eqRes;
    eigen_qp_lbA.block<20, 1>(6, 0) = neqRes_low;
    eigen_qp_ubA.block<6, 1>(0, 0) = eqRes;
    eigen_qp_ubA.block<20, 1>(6, 0) = neqRes_upp;

    Eigen::MatrixXd eigen_qp_H = Eigen::MatrixXd::Zero(QP_nv, QP_nv);
    Q2 = Eigen::MatrixXd::Identity(6, 6);
    Q1 = Eigen::MatrixXd::Identity(12, 12);
    // eigen_qp_H.block<6, 6>(0, 0) = Q2 * 2.0 * 1e7;
    // eigen_qp_H.block<12, 12>(6, 6) = Q1 * 2.0 * 1e1;
    eigen_qp_H.block<6, 6>(0, 0) = Q2 * 2.0 * 1e10;
    eigen_qp_H.block<12, 12>(6, 6) = Q1 * 2.0;

    // obj: (1/2)x'Hx+x'g
    // s.t. lbA<=Ax<=ubA
    //       lb<=x<=ub
//    qpOASES::real_t qp_H[QP_nv*QP_nv];
//    qpOASES::real_t qp_A[QP_nc*QP_nv];
//    qpOASES::real_t qp_g[QP_nv];
//    qpOASES::real_t qp_lbA[QP_nc];
//    qpOASES::real_t qp_ubA[QP_nc];
//    qpOASES::real_t xOpt_iniGuess[QP_nv];

    // std::cout<<"qp_A"<<std::endl;
    // std::cout<<eigen_qp_A_final<<std::endl;

    copy_Eigen_to_real_t(qp_H, eigen_qp_H, eigen_qp_H.rows(), eigen_qp_H.cols());
    copy_Eigen_to_real_t(qp_A, eigen_qp_A_final, eigen_qp_A_final.rows(), eigen_qp_A_final.cols());
    copy_Eigen_to_real_t(qp_lbA, eigen_qp_lbA, eigen_qp_lbA.rows(), eigen_qp_lbA.cols());
    copy_Eigen_to_real_t(qp_ubA, eigen_qp_ubA, eigen_qp_ubA.rows(), eigen_qp_ubA.cols());

    qpOASES::returnValue res;
    for (int i = 0; i < QP_nv; i++) {
        xOpt_iniGuess[i] = 0;
//        xOpt_iniGuess[i] =eigen_xOpt(i);
        qp_g[i] = 0;
    }
    nWSR = 200;
    cpu_time = timeStep;
//    QP_prob.reset();
    res = QP_prob.init(qp_H, qp_g, qp_A, NULL, NULL, qp_lbA, qp_ubA, nWSR, &cpu_time, xOpt_iniGuess);
    qpStatus = qpOASES::getSimpleStatus(res);
//    if (res==qpOASES::SUCCESSFUL_RETURN)
//        printf("WBC-QP: successful_return\n");
//    else if (res==qpOASES::RET_MAX_NWSR_REACHED)
//        printf("WBC-QP: max_nwsr\n");
//    else if (res==qpOASES::RET_INIT_FAILED)
//        printf("WBC-QP: init_failed\n");

    qpOASES::real_t xOpt[QP_nv];
    QP_prob.getPrimalSolution(xOpt);
    if (res == qpOASES::SUCCESSFUL_RETURN)
        for (int i = 0; i < QP_nv; i++)
            eigen_xOpt(i) = xOpt[i];

    eigen_ddq_Opt = ddq_final_kin;
    eigen_ddq_Opt.block<6, 1>(0, 0) += eigen_xOpt.block<6, 1>(0, 0);
    eigen_fr_Opt = Fr_ff + eigen_xOpt.block<12, 1>(6, 0);

    if (qpStatus != 0){
        Eigen::MatrixXd A_x;
        Eigen::VectorXd xOpt_iniGuess_m(QP_nv,1);
        for (int i=0; i < QP_nv;i++)
            xOpt_iniGuess_m(i) = xOpt_iniGuess[i];

    }

    Eigen::VectorXd tauRes;
    tauRes = dyn_M * eigen_ddq_Opt + dyn_Non - Jfe.transpose() * eigen_fr_Opt;

    tauJointRes = tauRes.block(6, 0, model_nv - 6, 1);
    // std::cout<<"qpRes_frOpt"<<std::endl;
    // std::cout<<eigen_fr_Opt.transpose()<<std::endl;
    std::cout<<"tauJointRes"<<std::endl;
    std::cout<<tauJointRes.transpose()<<std::endl;

    last_nWSR = nWSR;
    last_cpu_time = cpu_time;
}

void WBC_priority::computeDdq(Pin_KinDyn &pinKinDynIn) {
    // task definition
    /// -------- walk -------------
    {
        int id = kin_tasks.getId("static_Contact");
        kin_tasks.taskLib[id].errX = Eigen::VectorXd::Zero(6);
        kin_tasks.taskLib[id].derrX = Eigen::VectorXd::Zero(6);
        kin_tasks.taskLib[id].ddxDes = Eigen::VectorXd::Zero(6);
        kin_tasks.taskLib[id].dxDes = Eigen::VectorXd::Zero(6);
        kin_tasks.taskLib[id].kp = Eigen::MatrixXd::Identity(6, 6) * 0;
        kin_tasks.taskLib[id].kd = Eigen::MatrixXd::Identity(6, 6) * 0;
        kin_tasks.taskLib[id].J = Jc;
        kin_tasks.taskLib[id].dJ = dJc;
        kin_tasks.taskLib[id].W.diagonal() = Eigen::VectorXd::Ones(model_nv);

        // id = kin_tasks_walk.getId("RedundantJoints");
        // kin_tasks_walk.taskLib[id].errX = Eigen::VectorXd::Zero(5);
        // kin_tasks_walk.taskLib[id].errX(0) = 0 - q(21);
        // kin_tasks_walk.taskLib[id].errX(1) = 0 - q(22);
        // kin_tasks_walk.taskLib[id].errX(2) = 0 - q(23);
        // kin_tasks_walk.taskLib[id].errX(3) = 0 - q(24);
        // kin_tasks_walk.taskLib[id].errX(4) = 0 - q(25);
        // kin_tasks_walk.taskLib[id].derrX = Eigen::VectorXd::Zero(5);
        // kin_tasks_walk.taskLib[id].ddxDes = Eigen::VectorXd::Zero(5);
        // kin_tasks_walk.taskLib[id].dxDes = Eigen::VectorXd::Zero(5);
        // kin_tasks_walk.taskLib[id].kp = Eigen::MatrixXd::Identity(5, 5) * 200;
        // kin_tasks_walk.taskLib[id].kd = Eigen::MatrixXd::Identity(5, 5) * 20;
        // kin_tasks_walk.taskLib[id].J = Eigen::MatrixXd::Zero(5, model_nv);
        // kin_tasks_walk.taskLib[id].J(0, 20) = 1;
        // kin_tasks_walk.taskLib[id].J(1, 21) = 1;
        // kin_tasks_walk.taskLib[id].J(2, 22) = 1;
        // kin_tasks_walk.taskLib[id].J(3, 23) = 1;
        // kin_tasks_walk.taskLib[id].J(4, 24) = 1;
        // kin_tasks_walk.taskLib[id].dJ = Eigen::MatrixXd::Zero(5, model_nv);
        // kin_tasks_walk.taskLib[id].W.diagonal() = Eigen::VectorXd::Ones(model_nv);

        id = kin_tasks.getId("Roll_Pitch_Yaw");
        kin_tasks.taskLib[id].errX = Eigen::VectorXd::Zero(3);
        Eigen::Matrix3d desRot = eul2Rot(base_rpy_des(0), base_rpy_des(1), base_rpy_des(2));
        kin_tasks.taskLib[id].errX.block<3, 1>(0, 0) = diffRot(base_rot, desRot);
        kin_tasks.taskLib[id].derrX = Eigen::VectorXd::Zero(3);
        kin_tasks.taskLib[id].derrX.block<3, 1>(0, 0) = -dq.block<3, 1>(3, 0);
        kin_tasks.taskLib[id].ddxDes = Eigen::VectorXd::Zero(3);
        kin_tasks.taskLib[id].dxDes = Eigen::VectorXd::Zero(3);
        // kin_tasks.taskLib[id].kp = Eigen::MatrixXd::Identity(3, 3) * 2000;
        // kin_tasks.taskLib[id].kd = Eigen::MatrixXd::Identity(3, 3) * 100;
        kin_tasks.taskLib[id].kp = Eigen::MatrixXd::Identity(3, 3) * 100;
        kin_tasks.taskLib[id].kd = Eigen::MatrixXd::Identity(3, 3) * 25;
        kin_tasks.taskLib[id].J = J_base.block(3,0,3,model_nv);
        kin_tasks.taskLib[id].dJ = dJ_base.block(3,0,3,model_nv);
        kin_tasks.taskLib[id].W.diagonal() = Eigen::VectorXd::Ones(model_nv);

        id = kin_tasks.getId("PxPyPz");
        kin_tasks.taskLib[id].errX=Eigen::VectorXd::Zero(3);
        kin_tasks.taskLib[id].errX=base_pos_des - q.block(0,0,3,1);
        kin_tasks.taskLib[id].derrX = Eigen::VectorXd::Zero(3);
        kin_tasks.taskLib[id].ddxDes = Eigen::VectorXd::Zero(3);
        kin_tasks.taskLib[id].dxDes = Eigen::VectorXd::Zero(3);
        kin_tasks.taskLib[id].kp = Eigen::MatrixXd::Identity(3, 3) * 500; //100
        kin_tasks.taskLib[id].kd = Eigen::MatrixXd::Identity(3, 3) * 50;
        kin_tasks.taskLib[id].J = J_base.block(0,0,3,model_nv);
        kin_tasks.taskLib[id].dJ = dJ_base.block(0,0,3,model_nv);
        kin_tasks.taskLib[id].W.diagonal() = Eigen::VectorXd::Ones(model_nv);

        id = kin_tasks.getId("SwingLeg");
        kin_tasks.taskLib[id].errX = Eigen::VectorXd::Zero(6);
        kin_tasks.taskLib[id].errX.block<3, 1>(0, 0) = swing_F_fe_pos_des_W - F_fe_pos_sw_W;
        kin_tasks.taskLib[id].errX.block<3, 1>(3, 0) = swing_R_fe_pos_des_W - R_fe_pos_sw_W;
        std::cout<<"swing_F_fe_pos_des_W: "<<swing_F_fe_pos_des_W.transpose()<<"  "<<swing_R_fe_pos_des_W.transpose()<<std::endl;
        std::cout<<"F_fe_pos_sw_W: "<<F_fe_pos_sw_W.transpose()<<"  "<<R_fe_pos_sw_W.transpose()<<std::endl;   
        std::cout<<"errX: "<<kin_tasks.taskLib[id].errX.transpose()<<std::endl;
        kin_tasks.taskLib[id].derrX = Eigen::VectorXd::Zero(6);
//        kin_tasks.taskLib[id].derrX=-Jsw*dq;
        kin_tasks.taskLib[id].ddxDes = Eigen::VectorXd::Zero(6);
        kin_tasks.taskLib[id].dxDes = Eigen::VectorXd::Zero(6);
        kin_tasks.taskLib[id].kp = Eigen::MatrixXd::Identity(6, 6) * 2000;
        // kin_tasks.taskLib[id].kp.block<1, 1>(2, 2) = kin_tasks.taskLib[id].kp.block<1, 1>(2, 2) * 0.1;
        kin_tasks.taskLib[id].kd = Eigen::MatrixXd::Identity(6, 6) * 20;
        // kin_tasks.taskLib[id].kp = Eigen::MatrixXd::Identity(6, 6) * 10;
        // // kin_tasks.taskLib[id].kp.block<1, 1>(2, 2) = kin_tasks.taskLib[id].kp.block<1, 1>(2, 2) * 0.1;
        // kin_tasks.taskLib[id].kd = Eigen::MatrixXd::Identity(6, 6) * 5;
        
        kin_tasks.taskLib[id].J = Jsw;
        kin_tasks.taskLib[id].dJ = dJsw;
        kin_tasks.taskLib[id].W.diagonal() = Eigen::VectorXd::Ones(model_nv);

    }

    kin_tasks.computeAll(des_delta_q, des_dq, des_ddq, dyn_M, dyn_M_inv, dq);

    delta_q_final_kin = kin_tasks.out_delta_q;
    dq_final_kin = kin_tasks.out_dq;
    ddq_final_kin = kin_tasks.out_ddq;

    // final WBC output collection
}

void WBC_priority::copy_Eigen_to_real_t(qpOASES::real_t *target, const Eigen::MatrixXd &source, int nRows, int nCols) {
    int count = 0;

    for (int i = 0; i < nRows; i++) {
        for (int j = 0; j < nCols; j++) {
            target[count++] = isinf(source(i, j)) ? qpOASES::INFTY : source(i, j);
        }
    }
}

void WBC_priority::setQini(const Eigen::VectorXd &qIniDesIn, const Eigen::VectorXd &qIniCurIn) {
    qIniDes = qIniDesIn;
    qIniCur = qIniCurIn;
}