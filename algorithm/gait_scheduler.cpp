/*
This is part of OpenLoong Dynamics Control, an open project for the control of biped robot,
Copyright (C) 2024 Humanoid Robot (Shanghai) Co., Ltd, under Apache 2.0.
Feel free to use in any purpose, and cite OpenLoong-Dynamics-Control in any style, to contribute to the advancement of the community.
 <https://atomgit.com/openloong/openloong-dyn-control.git>
 <web@openloong.org.cn>
*/

#include "gait_scheduler.h"

// Note: no double-support here, swing time always equals to stance time
GaitScheduler::GaitScheduler(double tSwingIn, double dtIn) {
    tSwing = tSwingIn;
    dt = dtIn;
    phi = 0;
    isIni = false;
    legState = DataBus::FrSt;
}

void GaitScheduler::dataBusRead(const DataBus &robotState) {
    model_nv = robotState.model_nv;
    torJoint = Eigen::VectorXd::Zero(model_nv - 6);
    for (int i = 0; i < model_nv - 6; i++) {
        torJoint[i] = robotState.motors_tor_cur[i];
    }
    dyn_M = robotState.dyn_M;
    dyn_Non = robotState.dyn_Non;
    J_fr = robotState.J_fr;
    J_fl = robotState.J_fl;
    J_rr = robotState.J_rr;
    J_rl = robotState.J_rl;

    dJ_fr = robotState.dJ_fr;
    dJ_fl = robotState.dJ_fl;
    dJ_rr = robotState.dJ_rr;
    dJ_rl = robotState.dJ_rl;

    hip_fr_pos_W = robotState.hip_fr_pos_W;
    hip_fl_pos_W = robotState.hip_fl_pos_W;
    hip_rr_pos_W = robotState.hip_rr_pos_W;
    hip_rl_pos_W = robotState.hip_rl_pos_W;
    fe_fr_pos_W = robotState.fe_fr_pos_W;
    fe_fl_pos_W = robotState.fe_fl_pos_W;
    fe_rr_pos_W = robotState.fe_rr_pos_W;
    fe_rl_pos_W = robotState.fe_rl_pos_W;
    fe_l_rot_W = robotState.fe_l_rot_W;
    fe_r_rot_W = robotState.fe_r_rot_W;
    dq = robotState.dq;
}

void GaitScheduler::dataBusWrite(DataBus &robotState) {
    robotState.tSwing = tSwing;
    robotState.swingStartPos_F_W = swingStartPos_F_W;
    robotState.swingStartPos_R_W = swingStartPos_R_W;
    robotState.stanceDesPos_F_W = stanceStartPos_F_W;
    robotState.stanceDesPos_R_W = stanceStartPos_R_W;
    robotState.posHip_F_W = posHip_F_W;
    robotState.posHip_R_W = posHip_R_W;
    robotState.theta0 = theta0;
    robotState.legState = legState;
    robotState.legStateNext = legStateNext;
    robotState.phi = phi;
    robotState.FL_est = F_LFest;
    robotState.FR_est = F_RFest;

    if (legState == DataBus::FlSt) {
        robotState.stance_F_fe_pos_cur_W = fe_fl_pos_W;
        robotState.stance_R_fe_pos_cur_W = fe_rr_pos_W;
    } else {
        robotState.stance_F_fe_pos_cur_W = fe_fr_pos_W;
        robotState.stance_R_fe_pos_cur_W = fe_rl_pos_W;
    }
    std::cout<<"legState: "<<legState<<std::endl;
}

void GaitScheduler::step() {
    Eigen::VectorXd tauAll;
    tauAll = Eigen::VectorXd::Zero(model_nv);
    tauAll.block(6, 0, model_nv - 6, 1) = torJoint;
    F_LFest = -pseudoInv_SVD(J_fr * dyn_M.inverse() * J_fr.transpose()) *
            (J_fr * dyn_M.inverse() * (tauAll - dyn_Non) + dJ_fr * dq);//pseudoInv_SVD：奇异值分解法求伪逆；
    F_RFest = -pseudoInv_SVD(J_fl * dyn_M.inverse() * J_fl.transpose()) *
            (J_fl * dyn_M.inverse() * (tauAll - dyn_Non) + dJ_fl * dq);
    F_LHest = -pseudoInv_SVD(J_rr * dyn_M.inverse() * J_rr.transpose()) *
            (J_rr * dyn_M.inverse() * (tauAll - dyn_Non) + dJ_rr * dq);//pseudoInv_SVD：奇异值分解法求伪逆；
    F_RHest = -pseudoInv_SVD(J_rl * dyn_M.inverse() * J_rl.transpose()) *
            (J_rl * dyn_M.inverse() * (tauAll - dyn_Non) + dJ_rl * dq);

    double dPhi = 1.0 / tSwing * dt;

    phi += dPhi;
    if (!isIni) {
        isIni = true;
        if (legState == DataBus::FlSt) {
            swingStartPos_F_W = fe_fr_pos_W;
            swingStartPos_R_W = fe_rl_pos_W;
            stanceStartPos_F_W = fe_fl_pos_W;
            stanceStartPos_R_W = fe_rr_pos_W;
        } else {
            swingStartPos_F_W = fe_fl_pos_W;
            swingStartPos_R_W = fe_rr_pos_W;
            stanceStartPos_F_W = fe_fr_pos_W;
            stanceStartPos_R_W = fe_rl_pos_W;
        }
    }

    // if (legState == DataBus::FlSt && FRest[2] >= FzThrehold && phi>0.6) {
    //     legState = DataBus::RSt;
    //     swingStartPos_W = fe_l_pos_W;
    //     stanceStartPos_W = fe_r_pos_W;
    //     phi = 0;
    // } else if (legState == DataBus::RSt && FLest[2] >= FzThrehold && phi>0.6) {
    //     legState = DataBus::LSt;
    //     swingStartPos_W = fe_r_pos_W;
    //     stanceStartPos_W = fe_l_pos_W;
    //     phi = 0;
    // }

    if (legState == DataBus::FlSt && phi>=1.0) {
        legState = DataBus::FrSt;
        swingStartPos_F_W = fe_fl_pos_W;
        swingStartPos_R_W = fe_rr_pos_W;
        stanceStartPos_F_W = fe_fr_pos_W;
        stanceStartPos_R_W = fe_rl_pos_W;
        phi = 0;
    } else if (legState == DataBus::FrSt && phi>=1.0) {
        legState = DataBus::FlSt;
        swingStartPos_F_W = fe_fr_pos_W;
        swingStartPos_R_W = fe_rl_pos_W;
        stanceStartPos_F_W = fe_fl_pos_W;
        stanceStartPos_R_W = fe_rr_pos_W;
        phi = 0;
    }

    if (phi >= 1) {
        phi = 1;
    }
    if (legState == DataBus::FlSt) {
        posHip_F_W = hip_fr_pos_W;
        posHip_R_W = hip_rl_pos_W;
        theta0 = -3.1415 * 0.5;
        legStateNext = DataBus::FrSt;
    } else {
        posHip_F_W = hip_fl_pos_W;
        posHip_R_W = hip_rr_pos_W;
        theta0 = 3.1415 * 0.5;
        legStateNext = DataBus::FlSt;
    }
}



















