/*
This is part of OpenLoong Dynamics Control, an open project for the control of biped robot,
Copyright (C) 2024 Humanoid Robot (Shanghai) Co., Ltd, under Apache 2.0.
Feel free to use in any purpose, and cite OpenLoong-Dynamics-Control in any style, to contribute to the advancement of the community.
 <https://atomgit.com/openloong/openloong-dyn-control.git>
 <web@openloong.org.cn>
*/
#pragma once

#include "data_bus.h"
#include <Eigen/Dense>
#include "useful_math.h"

class GaitScheduler {
public:
    bool isIni{false};
    double phi{0};
    double tSwing{0.3};
    double dt{0.001};
    double FzThrehold{100};
    DataBus::LegState legState, legStateNext;
    GaitScheduler(double tSwingIn, double dtIn);
    void dataBusRead(const DataBus &robotState);
    void dataBusWrite(DataBus &robotState);
    void step();
    void stop();
    Eigen::VectorXd F_LFest, F_RFest, F_LHest,F_RHest;
    Eigen::VectorXd torJoint;
private:
    Eigen::VectorXd fe_fr_pos_W, fe_fl_pos_W, fe_rr_pos_W, fe_rl_pos_W, posHip_F_W, posHip_R_W, posST_W, dq;
    Eigen::VectorXd hip_fr_pos_W, hip_fl_pos_W, hip_rr_pos_W, hip_rl_pos_W;
    Eigen::VectorXd stanceStartPos_F_W, stanceStartPos_R_W, swingStartPos_F_W, swingStartPos_R_W;
    Eigen::MatrixXd fe_r_rot_W, fe_l_rot_W;
    Eigen::MatrixXd dyn_M, dyn_Non, J_fr, J_fl, J_rr, J_rl, dJ_fr, dJ_fl, dJ_rr, dJ_rl;
    double theta0;
    int model_nv;
};


