/*
This is part of OpenLoong Dynamics Control, an open project for the control of biped robot,
Copyright (C) 2024 Humanoid Robot (Shanghai) Co., Ltd, under Apache 2.0.
Feel free to use in any purpose, and cite OpenLoong-Dynamics-Control in any style, to contribute to the advancement of the community.
 <https://atomgit.com/openloong/openloong-dyn-control.git>
 <web@openloong.org.cn>
*/
#include "foot_placement.h"
#include "bezier_1D.h"

void FootPlacement::dataBusRead(DataBus &robotState) {
    posStart_F_W=robotState.swingStartPos_F_W;
    posStart_R_W=robotState.swingStartPos_R_W;
    desV_W=robotState.js_vel_des;
    desWz_W=robotState.js_omega_des(2);
    curV_W=robotState.base_rot*robotState.dq.block<3,1>(0,0);
    phi=robotState.phi;
    hipPos_F_W=robotState.posHip_F_W;
    hipPos_R_W=robotState.posHip_R_W;
    base_pos=robotState.base_pos;
    tSwing= robotState.tSwing;
    theta0=robotState.theta0;
    yawCur=robotState.rpy[2];
    omegaZ_W=robotState.base_omega_W(2);
    hip_width=robotState.width_hips;
    legState=robotState.legState;
}

void FootPlacement::dataBusWrite(DataBus &robotState) {
    //robotState.swing_fe_rpy_des_W<<0,0,robotState.base_rpy_des(2); // WARNING! ThetaZ!
    robotState.swing_F_fe_pos_des_W<<pDesCur_F[0],pDesCur_F[1],pDesCur_F[2];
    robotState.swing_R_fe_pos_des_W<<pDesCur_R[0],pDesCur_R[1],pDesCur_R[2];
}
void FootPlacement::getSwingPos() {
    Eigen::Matrix3d KP;
    KP.setZero();
    KP(0,0)=kp_vx;KP(1,1)=kp_vy;KP(2,2)=0;

    // for linear velocity
    posDes_F_W=hipPos_F_W+KP*(desV_W-curV_W)*(-1)+0.5*tSwing*curV_W+
            curV_W*(1-phi)*tSwing;

    posDes_R_W=hipPos_R_W+KP*(desV_W-curV_W)*(-1)+0.5*tSwing*curV_W+
            curV_W*(1-phi)*tSwing;

    // for angular veloctity
    // double thetaF;
    // thetaF=yawCur+theta0+omegaZ_W*(1-phi)*tSwing+0.5*omegaZ_W*tSwing+kp_wz*(omegaZ_W-desWz_W);
    // posDes_W(0)+=0.5*hip_width* (cos(thetaF)-cos(yawCur+theta0));
    // posDes_W(1)+=0.5*hip_width* (sin(thetaF)-sin(yawCur+theta0));

    //posDes_F_W(0)=posDes_F_W(0) - 0.02;
    // posDes_W(2)=STPos_W(2)-0.04;
    posDes_F_W(2)=base_pos(2)-legLength;
    posDes_R_W(2)=base_pos(2)-legLength;

    double yOff=0.085; // positive for moving the leg inside
    if (legState==DataBus::FlSt){
        posDes_F_W(1)-=yOff;
        posDes_R_W(1)+=yOff;
    }
    else if (legState==DataBus::FrSt){
        posDes_F_W(1)+=yOff;
        posDes_R_W(1)-=yOff;
    }

    // std::cout<<"hipPos_F_W: "<<hipPos_F_W.transpose()<<std::endl;
    // std::cout<<"hipPos_R_W: "<<hipPos_R_W.transpose()<<std::endl;
    // std::cout<<"posDes_F_W: "<<posDes_F_W.transpose()<<std::endl;
    // std::cout<<"posDes_R_W: "<<posDes_R_W.transpose()<<std::endl;

    // cycloid trajectories
    if (phi < 1.0){
        pDesCur_F[0]=posStart_F_W(0)+(posDes_F_W(0)-posStart_F_W(0))/(2*3.1415)*(2*3.1415*phi-sin(2*3.1415*phi));
        pDesCur_F[1]=posStart_F_W(1)+(posDes_F_W(1)-posStart_F_W(1))/(2*3.1415)*(2*3.1415*phi-sin(2*3.1415*phi));

        pDesCur_R[0]=posStart_R_W(0)+(posDes_R_W(0)-posStart_R_W(0))/(2*3.1415)*(2*3.1415*phi-sin(2*3.1415*phi));
        pDesCur_R[1]=posStart_R_W(1)+(posDes_R_W(1)-posStart_R_W(1))/(2*3.1415)*(2*3.1415*phi-sin(2*3.1415*phi));
    }
    // pDesCur[2]=posStart_W(2)+stepHeight*0.5*(1-cos(2*3.1415*phi))+(posDes_W(2)-posStart_W(2))/(2*3.1415)*(2*3.1415*phi-sin(2*3.1415*phi));
    pDesCur_F[2]=posStart_F_W(2)+Trajectory(tSwing/2, stepHeight, posDes_F_W(2)-posStart_F_W(2));
    pDesCur_R[2]=posStart_R_W(2)+Trajectory(tSwing/2, stepHeight, posDes_R_W(2)-posStart_R_W(2));

    //std::cout<<"pDesCur_F: "<<pDesCur_F[0]<<"  "<<pDesCur_F[1]<<"  "<<pDesCur_F[2]<<std::endl;

    for(int i=0; i<2; i++){
        pDesCur_F[i]=posStart_F_W(i);
        pDesCur_R[i]=posStart_R_W(i);
    }

}

double    FootPlacement::Trajectory(double phase, double des1, double des2){
    Bezier_1D Bswpid;
    double  para0 = 5, para1 = 3;
    for (int i = 0; i < para0; i ++)
        Bswpid.P.push_back(0.0);
    for (int i = 0; i < para1; i ++)
        Bswpid.P.push_back(1.0);

    double Bsw1, Bsw2;
    Bsw1 = Bswpid.getOut(phi/phase);
    Bsw2 = Bswpid.getOut((1.4 - phi)/(1.4 - phase));

    double output;
    if (phi < phase)
        output = des1*Bsw1;
    else{
        if (Bsw2 > 0)
            output = des1 * Bsw2 + des2 * (1.0 - Bsw2);
        else
            output = des2;
    }

    return output;
}

