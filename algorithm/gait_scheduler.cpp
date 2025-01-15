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
    hip_fr_pos_B = robotState.hip_fr_pos_B;
    hip_fl_pos_B = robotState.hip_fl_pos_B;
    hip_rr_pos_B = robotState.hip_rr_pos_B;
    hip_rl_pos_B = robotState.hip_rl_pos_B;
    fe_fr_pos_W = robotState.fe_fr_pos_W;
    fe_fl_pos_W = robotState.fe_fl_pos_W;
    fe_rr_pos_W = robotState.fe_rr_pos_W;
    fe_rl_pos_W = robotState.fe_rl_pos_W;
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
    robotState.posHip_F_B = posHip_F_B;
    robotState.posHip_R_B = posHip_R_B;
    robotState.F_theta0 = F_theta0;
    robotState.R_theta0 = R_theta0;
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
}

void GaitScheduler::step() {
    Eigen::VectorXd tauAll;
    tauAll = Eigen::VectorXd::Zero(model_nv);
    tauAll.block(6, 0, model_nv - 6, 1) = torJoint;
    F_LFest = -pseudoInv_SVD(J_fr.block(0,0,3,model_nv) * dyn_M.inverse() * J_fr.block(0,0,3,model_nv).transpose()) *
            (J_fr.block(0,0,3,model_nv) * dyn_M.inverse() * (tauAll - dyn_Non) + dJ_fr.block(0,0,3,model_nv) * dq);//pseudoInv_SVD：奇异值分解法求伪逆；
    F_RFest = -pseudoInv_SVD(J_fl.block(0,0,3,model_nv) * dyn_M.inverse() * J_fl.block(0,0,3,model_nv).transpose()) *
            (J_fl.block(0,0,3,model_nv) * dyn_M.inverse() * (tauAll - dyn_Non) + dJ_fl.block(0,0,3,model_nv) * dq);
    F_LHest = -pseudoInv_SVD(J_rr.block(0,0,3,model_nv) * dyn_M.inverse() * J_rr.block(0,0,3,model_nv).transpose()) *
            (J_rr.block(0,0,3,model_nv) * dyn_M.inverse() * (tauAll - dyn_Non) + dJ_rr.block(0,0,3,model_nv) * dq);//pseudoInv_SVD：奇异值分解法求伪逆；
    F_RHest = -pseudoInv_SVD(J_rl.block(0,0,3,model_nv) * dyn_M.inverse() * J_rl.block(0,0,3,model_nv).transpose()) *
            (J_rl.block(0,0,3,model_nv) * dyn_M.inverse() * (tauAll - dyn_Non) + dJ_rl.block(0,0,3,model_nv) * dq);

    double dPhi = 1.0 / tSwing * dt;

    // std::cout<<"Fest:"<<F_LFest.transpose()<<"  "<<F_RFest.transpose()<<"  "<<F_LHest.transpose()<<"  "<<F_RHest.transpose()<<std::endl;

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

    // if (legState == DataBus::FlSt && (F_RFest[2]>=FzThrehold || F_LHest[2]>=FzThrehold) && phi>=0.6) {
    if (legState == DataBus::FlSt && phi>=1) {
        legState = DataBus::FrSt;
        swingStartPos_F_W = fe_fl_pos_W;
        swingStartPos_R_W = fe_rr_pos_W;
        stanceStartPos_F_W = fe_fr_pos_W;
        stanceStartPos_R_W = fe_rl_pos_W;
        phi = 0;
    } 
    // else if (legState == DataBus::FrSt  && (F_RHest[2]>=FzThrehold || F_LFest[2]>=FzThrehold) && phi>=0.6) {
    else if (legState == DataBus::FrSt  && phi>=1) {
        legState = DataBus::FlSt;
        swingStartPos_F_W = fe_fr_pos_W;
        swingStartPos_R_W = fe_rl_pos_W;
        stanceStartPos_F_W = fe_fl_pos_W;
        stanceStartPos_R_W = fe_rr_pos_W;
        phi = 0;
    }

    std::cout<<"phi: "<<phi<<std::endl;

    if (phi >= 1) {
        phi = 1;
    }
    if (legState == DataBus::FlSt) {
        posHip_F_W = hip_fr_pos_W;
        posHip_R_W = hip_rl_pos_W;
        posHip_F_B = hip_fr_pos_B;
        posHip_R_B = hip_rl_pos_B;
        F_theta0 = -0.53;
        R_theta0 = 3.1415-0.53;
        legStateNext = DataBus::FrSt;
    } else {
        posHip_F_W = hip_fl_pos_W;
        posHip_R_W = hip_rr_pos_W;
        posHip_F_B = hip_fl_pos_B;
        posHip_R_B = hip_rr_pos_B;

        F_theta0 = 0.53;
        R_theta0 = -3.1415+0.53;
        legStateNext = DataBus::FlSt;
    }
}



















