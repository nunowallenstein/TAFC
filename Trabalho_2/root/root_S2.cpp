#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <cmath>
  using namespace std;
void root_S2()
{


  static const Double_t length=15;
  static const Int_t N=10;
  static const Int_t N_iter=2000;
  static const Double_t dt =0.005;

  //Vamos por numa file os tempos, as posições e as velocidades
  ofstream file;
  file.open("output.txt");
  

   
  // Parametros

  Double_t sig=1;
  Double_t m=4* M_PI;
  Double_t n_0=1;
  Double_t w_p=sqrt(4 * M_PI*pow(sig,2)/m);

  
  Double_t x_i=1;
  Double_t v_i=0;


  //ROOT

  TGraph *gr[N];
  TH1D *h1[1];
  h1[0]=new TH1D("h01","Velocities histogram",1000,-3,3);

 
  
  TCanvas *c1 = new TCanvas("c1","Positions", 750, 650);
  TCanvas *c2 = new TCanvas("c2","Velocities", 750, 650);

  c1->Divide(1,1);
  c2->Divide(1,1);


    //Declaração dos parametros
  Double_t *x_eq;
  x_eq=new Double_t [N];
  Double_t **x=new Double_t* [N];
  Double_t **v=new Double_t* [N];
  Double_t *t=new Double_t[N_iter];

  

  //Inicialização
  file << "0 "; 
  for (int i=0;i<N;i++)
    {  
      x[i]=new Double_t[N_iter];
      v[i]=new Double_t[N_iter];
      
      x_eq[i]=(length/(2*N))+ (i * length/N);
      x[i][0]=x_eq[i]+0.5*pow((+1),i);
      // x[i][0]=x_eq[i];      
      v[i][0]=v_i;
      file << x[i][0] << " " << v[i][0] << " ";
      if (x[i][0]>length || x[i][0]<0)
	{

	  cout <<"It's not possible to have initial values out of the bounds of the plasma, please make sure the initial positions are within the interval of [0,"<<length<<"]"<<endl;
	  return;
	}
    }

  
  file << endl;
  Double_t dx,dv;
  // v[N-1][0]=-2;

  for(int i=1;i<N_iter;i++)
    {
      Double_t taux=i*dt;
      t[i]=taux;
      file << taux << " ";
      /*
	for(int j=0;j<N;j++)
	{
	  
	dx_old=x[j][i-1]-x_eq[j][i-1];
	// cout << dx<<endl;
	dv_old=v[j][i-1]-v_eq[j][i-1];
	  
	//Força=-dx
	dv=dv_old-dx_old*dt;	  
	//	  incremento nas velocidades
	v_eq[j][i]=v_eq[j][i-1];
	v[j][i]=dv+v_eq[j][i];
	 	  

	//incremento nas posições
	dx=dx_old+dv*dt;
	 
	x_eq[j][i]=x_eq[j][i-1];
	x[j][i]=x_eq[j][i]+dx;
		  
	 
	//	 cout << x[j][i] <<endl;   
	file << x[j][i] << " " << v[j][i] << " ";	  
	}
      */

      
      for(int j=0;j<N;j++)
	{
	  dx=(+1)*(x[j][i-1]-x_eq[j]);
	  v[j][i]=v[j][i-1]*cos(w_p*dt)-w_p*dx*sin(w_p)*dt;	
	  x[j][i]=x[j][i-1]+v[j][i-1]*sin(w_p*dt)-dx*(1-cos(w_p*dt));

	
	  //Condições de Fronteira periódicas
	  if (x[j][i]<0)
	  {
	  x[j][i]=length;
	  x_eq[j]=x_eq[j] + length;	 
	  }

	  if (x[j][i]>length)
	  {
	  x[j][i]=0;
	  x_eq[j]=x_eq[j]-length;	 
	  }
	  file << x[j][i] << " " << v[j][i] << " ";
        
	}
      file << endl;
    }
  file.close();

  
  //Output

  c1->cd(1);
  for(int i=0;i<N;i++)
    {
      gr[i]=new TGraph(N_iter,x[i],t);

             auto *axis = gr[i]->GetXaxis();
	     gr[i]->SetMarkerStyle(1);
 
         axis->SetLimits(0-2,length+2);                 // along X
	 /*
	 gr[i]->GetHistogram()->SetMaximum(10);   // along          
      gr[i]->GetHistogram()->SetMinimum(-10);  //   Y     
	 */
      if (i==0)
	{
	  gr[i]->Draw("ap");
     
	}
      else
	gr[i]->Draw("p1 same");
    }


  //Integração trapézio
  Double_t *Int_V;
  Int_V=new Double_t[N];
  for (Int_t i=0;i<N;i++)
    {
      for(Int_t j=0;j<N_iter;j++)
	{
	  Int_V[i]+=v[i][j]*dt+(v[i][j+1]-v[i][j])*dt/2;
	}
      Int_V[i]=Int_V[i]/(N_iter*dt);
      h1[0]->Fill(Int_V[i]);
    }

   c2->cd(1);
   h1[0]->Draw();
  delete x_eq;
  delete[] x;
  delete[] v;

  
  
}
