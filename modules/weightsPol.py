import sys
import math
import ROOT
from array import array
from podio import root_io
import edm4hep 
 
def newAtau(TauP4, MesonP4,Type,New_Atau):
 
    if (Type!=0 and Type!=1 and Type!=10): # muons/electrons and others not implemented yet
         return 1 

    sin2theta_effective= 0.2312
    gv_ga=  -1 + 4 *sin2theta_effective
    Ae_sm=  2* gv_ga / (1+gv_ga*gv_ga)
    Atau_sm= Ae_sm    

    Ae= Ae_sm 

    # Polarization depends on cos(Theta):
    costheta_tau=math.cos(TauP4.Theta()) # this is the theta of the Tau, not the meson 
    Ptau_sm= - (Atau_sm * (1+  costheta_tau*costheta_tau) + 2*Ae_sm*costheta_tau) / (1+costheta_tau*costheta_tau + 2*Ae_sm*Atau_sm*costheta_tau)
    Pnew = - ( New_Atau   * (1+  costheta_tau*costheta_tau) + 2*Ae*costheta_tau) / (1+costheta_tau*costheta_tau + 2*Ae* New_Atau *costheta_tau)

    alpha= 1 
    if (Type==0):
        alpha= 1
    elif (Type==1):
        alpha= 0.46
    elif (Type==10):
        alpha= 0.12

    # rest frame of the tau?

    boost=ROOT.TVector3()
    boost=TauP4.BoostVector()

    meson_TauRes=ROOT.TLorentzVector()
    meson_TauRes.SetPxPyPzE(MesonP4.Px(),MesonP4.Py(),MesonP4.Pz(),MesonP4.E())
    meson_TauRes.Boost(-boost)

    # theta angle between meson and tau 
    # θ = cos-1 [ (a · b) / (|a| |b|) ].
    v1 = meson_TauRes.Vect()
    v2 = TauP4.Vect()
    theta_meson= v1.Angle(v2)
    z = math.cos( theta_meson)
    #print (theta_meson,z)

    weight_Pnew = (1+ alpha*Pnew*z) /(1+alpha*Ptau_sm*z)

    return weight_Pnew



def newAtauRHO(TauP4, RhoP4,beamE, TauConst, Type,New_Atau,sin2theta_effective= 0.2312):
 
    if (Type!=1): # this is for RHOs
         weight=newAtau(TauP4, RhoP4,Type,New_Atau)                  
         return weight  
 
    #sin2theta_effective= 0.2312
    gv_ga=  -1 + 4 *sin2theta_effective
    Ae_sm=  2* gv_ga / (1+gv_ga*gv_ga)
    Atau_sm= Ae_sm    

    Ae= Ae_sm

    # Polarization depends on cos(Theta):
    costheta_tau=math.cos(TauP4.Theta()) # this is the theta of the Tau, not the meson 
    Ptau_sm= - (Atau_sm * (1+  costheta_tau*costheta_tau) + 2*Ae_sm*costheta_tau) / (1+costheta_tau*costheta_tau + 2*Ae_sm*Atau_sm*costheta_tau)
    Pnew = - ( New_Atau   * (1+  costheta_tau*costheta_tau) + 2*Ae*costheta_tau) / (1+costheta_tau*costheta_tau + 2*Ae* New_Atau *costheta_tau)

    alpha= 0.46

    # pion?

    # two constituents in this case
    #print (Type,len(TauConst),TauConst[0].getPDG())
    pion=TauConst[0]
    pi0=TauConst[1]
    if abs(pion.getPDG())!=211:
      pion=TauConst[1]
      pi0=TauConst[0]
    #print ("pion ",pion.getPDG()," pi0 ",pi0.getPDG())

    pionP4=ROOT.TLorentzVector()
    pionP4.SetXYZM(pion.getMomentum().x,pion.getMomentum().y,pion.getMomentum().z,pion.getMass())

    mtau=TauP4.M()
    mRho=RhoP4.M()

    # rest frame of the tau?

    boost=ROOT.TVector3()
    boost=TauP4.BoostVector()

    Rho_TauRes=ROOT.TLorentzVector()
    Rho_TauRes.SetPxPyPzE(RhoP4.Px(),RhoP4.Py(),RhoP4.Pz(),RhoP4.E())
    Rho_TauRes.Boost(-boost)

    boostRho=ROOT.TVector3()
    boostRho=RhoP4.BoostVector()

    Pion_RhoRes=ROOT.TLorentzVector()
    Pion_RhoRes.SetPxPyPzE(pionP4.Px(),pionP4.Py(),pionP4.Pz(),pionP4.E())
    Pion_RhoRes.Boost(-boostRho)

#    We have θ = cos-1 [ (a · b) / (|a| |b|) ].e.  
    v1 = Rho_TauRes.Vect()
    v2 = TauP4.Vect()
    theta_Rho= v1.Angle(v2)
    z = math.cos( theta_Rho) 
    #print (theta_Rho,z)
    
    v1 = Pion_RhoRes.Vect()
    v2 = RhoP4.Vect()
    beta= v1.Angle(v2)
    cosBeta = math.cos( beta)
    #print (theta_Rho,z)

    x=RhoP4.E()/TauP4.E()

    sqrts=beamE*2
 
    cosPsi= (x * (mtau*mtau + mRho*mRho) - 2*mRho*mRho)/((mtau*mtau-mRho*mRho)*math.sqrt(x*x-4*mRho*mRho/sqrts/sqrts))

    if cosPsi>1:
    #    print ('What happened?', cosPsi)
       cosPsi=1
    if cosPsi<-1:
    #    print ('What happened?', cosPsi)
       cosPsi=-1

    anglePsi=math.acos(cosPsi)


    # more complex version that takes into account the dependence on Q, B, Theta: 

     
    ratioMass2= mtau*mtau/mRho/mRho

    term_a= 2/3*((1-Ptau_sm*z)- ratioMass2*(1+Ptau_sm*z)) + ratioMass2* (1+Ptau_sm*z)
    term_b= -2/3*((1-Ptau_sm*z-ratioMass2*(1+Ptau_sm*z))*(3*cosPsi*cosPsi-1)/2-3/2*math.sqrt(ratioMass2)*Ptau_sm*math.sin(2*anglePsi) *math.sin(theta_Rho))*(3*cosBeta*cosBeta-1)/2
    den = term_a+term_b

    term_a= 2/3*((1-Pnew*z)- ratioMass2*(1+Pnew*z)) + ratioMass2* (1+Pnew*z)
    term_b= -2/3*((1-Pnew*z-ratioMass2*(1+Pnew*z))*(3*cosPsi*cosPsi-1)/2-3/2*math.sqrt(ratioMass2)*Pnew*math.sin(2*anglePsi) *math.sin(theta_Rho))*(3*cosBeta*cosBeta-1)/2

    num=term_a+term_b

    weight_Pnew =num/den

    return weight_Pnew




def newAtauRHO2(TauP4, RhoP4,pionP4,beamE, Type,New_Atau,sin2theta_effective= 0.2312):
 
    if (Type!=1): # this is for RHOs
         weight=newAtau(TauP4, RhoP4,Type,New_Atau)                  
         return weight  
 
    #sin2theta_effective= 0.2312
    gv_ga=  -1 + 4 *sin2theta_effective
    Ae_sm=  2* gv_ga / (1+gv_ga*gv_ga)
    Atau_sm= Ae_sm    

    Ae=  Ae_sm

    # Polarization depends on cos(Theta):
    costheta_tau=math.cos(TauP4.Theta()) # this is the theta of the Tau, not the meson 
    Ptau_sm= - (Atau_sm * (1+  costheta_tau*costheta_tau) + 2*Ae_sm*costheta_tau) / (1+costheta_tau*costheta_tau + 2*Ae_sm*Atau_sm*costheta_tau)
    Pnew = - ( New_Atau   * (1+  costheta_tau*costheta_tau) + 2*Ae*costheta_tau) / (1+costheta_tau*costheta_tau + 2*Ae* New_Atau *costheta_tau)

    # pion?

    mtau=TauP4.M()
    mRho=RhoP4.M()

    # rest frame of the tau?

    boost=ROOT.TVector3()
    boost=TauP4.BoostVector()

    Rho_TauRes=ROOT.TLorentzVector()
    Rho_TauRes.SetPxPyPzE(RhoP4.Px(),RhoP4.Py(),RhoP4.Pz(),RhoP4.E())
    Rho_TauRes.Boost(-boost)

    boostRho=ROOT.TVector3()
    boostRho=RhoP4.BoostVector()

    Pion_RhoRes=ROOT.TLorentzVector()
    Pion_RhoRes.SetPxPyPzE(pionP4.Px(),pionP4.Py(),pionP4.Pz(),pionP4.E())
    Pion_RhoRes.Boost(-boostRho)

#    We have θ = cos-1 [ (a · b) / (|a| |b|) ].e.  
    v1 = Rho_TauRes.Vect()
    v2 = TauP4.Vect()
    theta_Rho= v1.Angle(v2)
    z = math.cos( theta_Rho) 
    #print (theta_Rho,z)
    
    v1 = Pion_RhoRes.Vect()
    v2 = RhoP4.Vect()
    beta= v1.Angle(v2)
    cosBeta = math.cos( beta)
    #print (theta_Rho,z)

    x=RhoP4.E()/TauP4.E()

    sqrts=beamE*2
 
    cosPsi= (x * (mtau*mtau + mRho*mRho) - 2*mRho*mRho)/((mtau*mtau-mRho*mRho)*math.sqrt(x*x-4*mRho*mRho/sqrts/sqrts))

    if cosPsi>1:
       print ('What happened?', cosPsi)
       cosPsi=1
    if cosPsi<-1:
       print ('What happened?', cosPsi)
       cosPsi=-1

    anglePsi=math.acos(cosPsi)


    # more complex version that takes into account the dependence on Q, B, Theta: 

     
    ratioMass2= mtau*mtau/mRho/mRho

    term_a= 2/3*((1-Ptau_sm*z)- ratioMass2*(1+Ptau_sm*z)) + ratioMass2* (1+Ptau_sm*z)
    term_b= -2/3*((1-Ptau_sm*z-ratioMass2*(1+Ptau_sm*z))*(3*cosPsi*cosPsi-1)/2-3/2*math.sqrt(ratioMass2)*Ptau_sm*math.sin(2*anglePsi) *math.sin(theta_Rho))*(3*cosBeta*cosBeta-1)/2
    den = term_a+term_b

    term_a= 2/3*((1-Pnew*z)- ratioMass2*(1+Pnew*z)) + ratioMass2* (1+Pnew*z)
    term_b= -2/3*((1-Pnew*z-ratioMass2*(1+Pnew*z))*(3*cosPsi*cosPsi-1)/2-3/2*math.sqrt(ratioMass2)*Pnew*math.sin(2*anglePsi) *math.sin(theta_Rho))*(3*cosBeta*cosBeta-1)/2

    num=term_a+term_b

    weight_Pnew =num/den

    return weight_Pnew


