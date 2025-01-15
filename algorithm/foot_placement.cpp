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
    curV_W=robotState.dq.block<3,1>(0,0);
    phi=robotState.phi;

    hipPos_F_W=robotState.posHip_F_W;
    hipPos_R_W=robotState.posHip_R_W;
    hipPos_F_B=robotState.posHip_F_B;
    hipPos_R_B=robotState.posHip_R_B;

    // STPos_W=robotState.posST_W;
    base_pos=robotState.base_pos;
    tSwing= robotState.tSwing;

    F_theta0=robotState.F_theta0;
    R_theta0=robotState.R_theta0;
    yawCur=robotState.rpy[2];
    omegaZ_W=robotState.base_omega_W(2);
    hip_width=robotState.width_hips;
    legState=robotState.legState;
}

void FootPlacement::dataBusWrite(DataBus &robotState) {
    // robotState.swing_fe_rpy_des_W<<0,0,robotState.base_rpy_des(2); // WARNING! ThetaZ!
    robotState.swing_F_fe_pos_des_W<<pDesCur_F[0],pDesCur_F[1],pDesCur_F[2];
    robotState.swing_R_fe_pos_des_W<<pDesCur_R[0],pDesCur_R[1],pDesCur_R[2];
}
void FootPlacement::getSwingPos() {
    Eigen::Matrix<double,4,1> b;
    b.setZero();
    Eigen::Matrix<double,1,4> xNow;
    xNow<<1,phi,pow(phi,2),pow(phi,3);

    Eigen::Matrix3d KP, Rz;
    KP.setZero();
    KP(0,0)=kp_vx;KP(1,1)=kp_vy;KP(2,2)=0;
    Rz<<cos(yawCur),-sin(yawCur),0,
            sin(yawCur),cos(yawCur),0,
            0,0,1;
    KP=Rz*KP*Rz.transpose();

    bodyRadius = sqrt(pow(hipPos_F_B(0), 2) + pow(hipPos_F_B(1), 2));
    bodyRadius = 0.33775;

    // for linear velocity
    posDes_F_W=base_pos+KP*(desV_W-curV_W)*(-1)+0.5*tSwing*curV_W+
            curV_W*(1-phi)*tSwing;

    posDes_R_W=base_pos+KP*(desV_W-curV_W)*(-1)+0.5*tSwing*curV_W+
            curV_W*(1-phi)*tSwing;

    // for angular veloctity
    double F_thetaF;
    double R_thetaF;
    double thetaF;

    F_thetaF = yawCur+F_theta0+omegaZ_W*(1-phi)*tSwing+0.5*omegaZ_W*tSwing+kp_wz*(omegaZ_W-desWz_W);
    R_thetaF = yawCur+R_theta0+omegaZ_W*(1-phi)*tSwing+0.5*omegaZ_W*tSwing+kp_wz*(omegaZ_W-desWz_W);

    posDes_F_W(0)+=bodyRadius* cos(F_thetaF);
    posDes_F_W(1)+=bodyRadius* sin(F_thetaF);

    posDes_R_W(0)+=bodyRadius* cos(R_thetaF);
    posDes_R_W(1)+=bodyRadius* sin(R_thetaF);

    posDes_F_W(2)=0;
    posDes_R_W(2)=0;

    std::cout<<"posDes_F_W: "<<posDes_F_W.transpose()<<std::endl;
    std::cout<<"posDes_R_W: "<<posDes_R_W.transpose()<<std::endl;
    std::cout<<"desV_W: "<<desV_W<<std::endl;
    std::cout<<"legState: "<<legState<<std::endl;
    std::cout<<"yawCur: "<<yawCur<<std::endl;
    std::cout<<"bodyRadius: "<<bodyRadius<<std::endl;

    double xOff=-0.01; //-0.01; // foot-end position offset in x direction in body frame
    double yOff=0.01; // 0.01; // foot-end position offset in y direction in body frame, positive for moving the leg inside
    double zOff=-0.035; // foot-end position offset in z direction in world frame


    double xOff_W_F(0), xOff_W_R(0), yOff_W_F(0), yOff_W_R(0);
    if (legState==DataBus::FlSt) {
        xOff_W_F = cos(yawCur) * xOff - sin(yawCur) * yOff;
        yOff_W_F = sin(yawCur) * xOff + cos(yawCur) * yOff;
        xOff_W_R = cos(yawCur) * (-xOff) - sin(yawCur) * (-yOff+0.002);
        yOff_W_R = sin(yawCur) * (-xOff) + cos(yawCur) * (-yOff+0.002);
    } else if (legState==DataBus::FrSt){
        xOff_W_F = cos(yawCur) * xOff - sin(yawCur) * (-yOff);
        yOff_W_F = sin(yawCur) * xOff + cos(yawCur) * (-yOff);
        xOff_W_R = cos(yawCur) * (-xOff) - sin(yawCur) * (yOff-0.002);
        yOff_W_R = sin(yawCur) * (-xOff) + cos(yawCur) * (yOff-0.002);
    }

    // posDes_F_W(0)+= xOff_W_F;
    // posDes_F_W(1)+= yOff_W_F;
    // posDes_R_W(0)-= xOff_W_R;
    // posDes_R_W(1)-= yOff_W_R;

    // cycloid trajectories
    if (phi < 1.0){
        pDesCur_F[0]=posStart_F_W(0)+(posDes_F_W(0)-posStart_F_W(0))/(2*3.1415)*(2*3.1415*phi-sin(2*3.1415*phi));
        pDesCur_F[1]=posStart_F_W(1)+(posDes_F_W(1)-posStart_F_W(1))/(2*3.1415)*(2*3.1415*phi-sin(2*3.1415*phi));

        pDesCur_R[0]=posStart_R_W(0)+(posDes_R_W(0)-posStart_R_W(0))/(2*3.1415)*(2*3.1415*phi-sin(2*3.1415*phi));
        pDesCur_R[1]=posStart_R_W(1)+(posDes_R_W(1)-posStart_R_W(1))/(2*3.1415)*(2*3.1415*phi-sin(2*3.1415*phi));
    }
//    pDesCur[2]=posStart_W(2)+stepHeight*0.5*(1-cos(2*3.1415*phi))+(posDes_W(2)-posStart_W(2))/(2*3.1415)*(2*3.1415*phi-sin(2*3.1415*phi));

    pDesCur_F[2]=posStart_F_W(2)+Trajectory(tSwing/2, stepHeight, posDes_F_W(2)-posStart_F_W(2));
    pDesCur_R[2]=posStart_R_W(2)+Trajectory(tSwing/2, stepHeight, posDes_R_W(2)-posStart_R_W(2));
}

double FootPlacement::Trajectory(double phase, double hei, double len){
    Bezier_1D Bswpid;
    double para0=5, para1=3;
    for(int i=0; i<para0; i++){Bswpid.P.push_back(0.0);}
    for(int i=0; i<para1; i++){Bswpid.P.push_back(1.0);}

    double output;
    if(phi<phase){
        output=hei*Bswpid.getOut(phi/phase);
    }else{
        double s=Bswpid.getOut((1.4-phi)/(1.4-phase));
        if(s>0){
            output=hei*s +len*(1.0-s);
        }else{
            output=len;
        }
    }
    return output;
}

