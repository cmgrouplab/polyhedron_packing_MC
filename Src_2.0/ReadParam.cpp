
using namespace std;

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include "const.h" //for pi
#include "ReadParam.h"


ReadParam::ReadParam()
{
  cout <<" Please specify all the parameters for the packing program: "<<endl;

  cout <<" The number of stages -  Nstage := "; cin >> Nstage; cout<<endl;
  cout <<" The number of cycles per stage - Ncycle := "; cin >> Ncycle; cout<<endl;
  cout <<" The number of deformation per stage - Ndeform := "; cin>>Ndeform; cout<<endl;

  cout <<" The translation magnitude - TransMod := "; cin>>TransMod; cout<<endl;
  cout <<" The rotation magnitude (in unit pi) - RotMod := "; cin>>RotMod; RotMod = RotMod*pi; cout<<endl;
  cout <<" Normal strain magnitude - Global_StrainMod := "; cin>>Global_StrainMod; cout<<endl;
  cout <<" Shear strain magnitude - Shear_StrainMod := "; cin>>Shear_StrainMod; cout<<endl;

  cout <<" Strain rescale magnitude - Strain_Rescale := "; cin>>Strain_Rescale; cout<<endl;
  cout <<" Translation rescale magnitude - Trans_Rescale := "; cin>>Trans_Rescale; cout<<endl;
  cout <<" Rotation rescale magnitude - Rot_Rescale := "; cin>>Rot_Rescale; cout<<endl;

  cout <<" Probability of uphill moves - p_uphill := "; cin>>p_uphill; cout<<endl;
  cout <<" Probability of translation for a trail move - p_trans := "; cin>>p_trans; cout<<endl;

  cout <<" Rescale paramter for initial configuration - relax_ratio :="; cin>>relax_ratio; cout<<endl;
  cout <<" Starting density - Starting_Density := "; cin>>Starting_Density; cout<<endl;
  cout <<" Termination density - Terminate_Density := "; cin>>Terminate_Density; cout<<endl;

  cout<<" Rescale the packing or Exit if Iconfig.txt containing overlapping particles? Rescale - 0; Exit - 1"<<endl;
  cin >> check_option; //this also read in from the input.txt file
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  
  cout <<" The number of trial moves between each property collection - Npc := "; cin>>Npc; cout<<endl;

  cout <<" Collecting Pressure? 1 - Yes; 0 - No." <<endl;
  cin >> flag_comput_pressure; cout<<endl;
  if(flag_comput_pressure == 1)
    {
      cout <<" Pressure will be collected!" <<endl;
    }

  cout <<" Computing g2? 1 - Yes; 0 - No." <<endl;
  cin >> flag_comput_g2; cout<<endl;
  if(flag_comput_g2 == 1)
    {
      cout <<" g2 will be computed!" <<endl;
    }

  cout <<" Computing orientational order metric? 1 - Yes; 0 - No." <<endl;
  cin >> flag_comput_ordermetric; cout<<endl;
  if(flag_comput_ordermetric == 1)
    {
      cout <<" Order metric will be computed!" <<endl;
    }


  //cout <<" Near-Neighor-List (NNL) cut off distance ratio of d_cvmax - nnl_cut_dist :="; cin>>nnl_cut_dist; cout<<endl;
  //cout <<" Near-Neighor-List (NNL) threshold density (starting using NNL) - rho_nnl_th :="; cin>>rho_nnl_th; cout<<endl;
}

void ReadParam::Rescale_TransMod()
{
  TransMod = TransMod*Trans_Rescale;
}

void ReadParam::Rescale_RotMod()
{
  RotMod = RotMod*Rot_Rescale;
}

void ReadParam::Rescale_Strain()
{
  Global_StrainMod = Global_StrainMod*Strain_Rescale; //in the meantime, scale the strain magnitude
  Shear_StrainMod = Shear_StrainMod*Strain_Rescale;
}

int ReadParam::Get_Nstage()
{
  return Nstage;
}

int ReadParam::Get_Ncycle()
{
  return Ncycle;
}

int ReadParam::Get_Ndeform()
{
  return Ndeform;
}

double ReadParam::Get_TransMod()
{
  return TransMod;
}

double ReadParam::Get_RotMod()
{
  return RotMod;
}

double ReadParam::Get_Global_StrainMod()
{
  return Global_StrainMod;
}

double ReadParam::Get_Shear_StrainMod()
{
  return Shear_StrainMod;
}

double ReadParam::Get_p_uphill()
{
  return p_uphill;
}

double ReadParam::Get_p_trans()
{
  return p_trans;
}

double ReadParam::Get_relax_ratio()
{
  return relax_ratio;
}

double ReadParam::Get_Starting_Density()
{
  return Starting_Density;;
}

double ReadParam::Get_Terminate_Density()
{
  return Terminate_Density;;
}


int ReadParam::Get_Npc()
{
  return Npc;
}


int ReadParam::flag_pressure()
{
  return flag_comput_pressure;
}

int ReadParam::flag_g2()
{
  return flag_comput_g2;
}


int ReadParam::flag_ordermetric()
{
  return flag_comput_ordermetric;
}

int ReadParam::Get_check_option()
{
  return check_option;
}

double ReadParam::Get_Strain_Rescale()
{
  return Strain_Rescale;
}
