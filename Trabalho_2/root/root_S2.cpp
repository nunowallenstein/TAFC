
void root_S2()
{
 
  static const Double_t length=400;
  static const Int_t N=400;
  static const Int_t N_iter=6000;
  static const Double_t dt =0.005;



  //Vamos por numa file os tempos, as posições e as velocidades
  ofstream file;
  file.open("output.txt");


   
  // Parametros Físicos
  //parametros interessantes, m=1,e=1,w_p=3.530,v_th=1

  
  //Constantes
  Double_t m=1;
  Double_t e=1;
 

  //Variáveis Independentes
  Double_t w_p=3; //
  Double_t v_th=0.01;


  
  //Variáveis Dependentes
  Double_t kbT=(pow(v_th,2)*m)/2;
  Double_t n_0=m*pow(w_p,2)/(4*M_PI*pow(e,2));

  


  //Output Electric Field
  Int_t parts=300;     //partitions per sheet
      


  
  //ROOT
  //------------------------------------------------------------------------------------------------------------
  
  //Gráficos
  TGraph *gr[N];
  TGraph *E_gr[10];
  TGraph *Energy_gr;
  //Histograma
  TH1D *h1[1];
  h1[0]=new TH1D("h_1","Velocities histogram",10000,-10,10);

 
  //Canvas
  TCanvas *c1 = new TCanvas("c_pos","Positions", 750, 650);
  TCanvas *c2 = new TCanvas("c_Vel","Velocities", 750, 650);
  TCanvas *c3 = new TCanvas("c_E","Electric Field", 1000, 1000);
   TCanvas *c3_1 = new TCanvas("c_E_1","Electric Field", 1000, 1000);
  TCanvas *c4= new TCanvas("c_En","Energy", 750, 650);
  c1->Divide(1,1);
  c2->Divide(1,1);
  c3->Divide(4,1);
  c3_1->Divide(4,1);
  c4->Divide(1,1);
  //--------------------------------------------------------------------------------------------------------------
  //Declaração da dinâmica

  //Pos Equilibrio
  Double_t *x_eq;
  x_eq=new Double_t[N];

  //Posição,array 2d para a posição e tempo
  Double_t **x=new Double_t *[N];
  
  //Velocidade array 2d para a posição e tempo
  Double_t **v=new Double_t *[N];
 


  //----------------------------------------------------------------------------------------------------------------
  //Gaussiana

  
  Double_t *v_gauss;
  v_gauss=new Double_t[N];

  srand(time(NULL));
  
  //Box muller method to a non-unitary variance gaussian distribution http://alpheratz.net/Maple/GaussianDistribution/GaussianDistribution.pdf
   for(Int_t i=0;i<N;i++)
   {
    Double_t aux=((Double_t)rand()/(RAND_MAX));
   v_gauss[i]=v_th*sqrt(-2*log(aux))*cos(2* M_PI*aux);
  }

  /*  TRandom Gauss;
  for(Int_t i=0;i<N;i++)
    {
       v_gauss[i]=v_th*Gauss.Gaus(0,v_th);
      
    }
  */
  
  //-----------------------------------------------------------------------------------------------------------------
  //Inicialização
  file << "0 "; 
  for (int i=0;i<N;i++)
    {  

      x[i]=new Double_t[2];
      v[i]=new Double_t[2];
      
      x_eq[i]=(length/(2*N))+ (i * length/N);
      
      // x[i][0]=x_eq[i]+0.5*pow((+1),i);
      
      x[i][0]=x_eq[i];
    
      v[i][0]=v_gauss[i];
      if(i==N/2)
      	v[i][0]=10;
	 
	file << x[i][0]<< " " << v[i][0] << " ";
      //Caso as posições iniciais não estejam na caixa
      if (x[i][0]>length || x[i][0]<0)
	{
	  cout <<"It's not possible to have initial values out of the bounds of the plasma, please make sure the initial positions are within the interval of [0,"<<length<<"]"<<endl;
	  return;
	}
      
    }

  
  file << endl;
  //-----------------------------------------------------------------------------------------------------------------
  //Evolução no tempo
  
  
  Double_t dx,dv;
  Double_t dx_old,dv_old;
  Double_t dx_old_i,dx_old_j;
  Double_t v_old,x_eq_old,x_old;
  Double_t x_old_i,x_old_j;
  Double_t v_old_i,v_old_j;
  Double_t x_btw_i,x_btw_j;
  Double_t v_btw_i,v_btw_j;
    
  Double_t dtc1,dtc2,dt_f;
  
  //Variável trigger controla consoante a folha debaixo esteja em cima ou a decima em baixo


  //Crossings
  Int_t trigger_1=0;
  Int_t trigger_2=0;

  

  for(int i=1;i<N_iter;i++)
    {
      Double_t taux=i*dt;
      file << taux << " ";
          
      for(int j=0;j<N;j++)
	{
	  //guardar as variáveis x e y no tempo anterior, para aplicar a condição de crossing para j e j-1
	  x[j][1]=x[j][0];
	  v[j][1]=v[j][0];
	 
	  dx_old=(+1)*(x[j][1]-x_eq[j]);
	  if(j>0)
	    {
	      
	      x_old_i=x[j-1][1];
	      dx_old_i=x_old_i-x_eq[j-1];									    
	      x_old_j=x[j][1];
	      dx_old_j=x_old_j-x_eq[j];
	      v_old_i=v[j-1][1];
	      v_old_j=v[j][1];
	     
	    	  
	    }
	  
	  v[j][0]=v[j][0]*cos(w_p*dt)-w_p*dx_old*sin(w_p*dt);	
	  x[j][0]=x[j][0]+v[j][0]*sin(w_p*dt)-dx_old*(1-cos(w_p*dt));
	  
	  //Crossings
	  if(j>0)
	    {
	      if ((x[j-1][0]>x[j][0]))
		{
		  
		  //Avança-se dtc1
	
		  dtc1=dt*(x_old_j-x_old_i)/(x_old_j-x_old_i+x[j-1][0]-x[j][0]);
		 
		  
		  v_btw_i=v_old_i*cos(w_p*dtc1)-w_p*dx_old_i*sin(w_p*dtc1);	
		  v_btw_j=v_old_j*cos(w_p*dtc1)-w_p*dx_old_j*sin(w_p*dtc1);

		  x_btw_i=x_old_i+v_old_i*sin(w_p*dt)-dx_old_i*(1-cos(w_p*dtc1));
		  x_btw_j=x_old_j+v_old_j*sin(w_p*dt)-dx_old_j*(1-cos(w_p*dtc1));

		  
		  //Avança-se dtc2

		  dx_old_j=x_btw_j-x_eq[j];
		  dx_old_i=x_btw_i-x_eq[j-1];
		  
		  dtc2=dtc1*(x_old_j-x_old_i)/(x_old_j-x_old_i+x_btw_i-x_btw_j);

		 
		  v_old_i=v_btw_i*cos(w_p*dtc2)-w_p*dx_old_i*sin(w_p*dtc2);	
		  v_old_j=v_btw_j*cos(w_p*dtc2)-w_p*dx_old_j*sin(w_p*dtc2);

		  x_old_i=x_btw_i+v_btw_i*sin(w_p*dtc2)-dx_old_i*(1-cos(w_p*dtc2));
		  x_old_j=x_btw_j+v_btw_j*sin(w_p*dtc2)-dx_old_j*(1-cos(w_p*dtc2));

		  //Troca-se as velocidades
		  v_old=v_old_i;
		  v_old_i=v_old_j;
		  v_old_j=v_old;
		 
		  //Troca-se os centros de equilibrio
		  x_eq_old=x_eq[j-1];
		  x_eq[j-1]=x_eq[j];
		  x_eq[j]=x_eq_old;

		  //Propagation for the remainder of the timestep
		  dt_f=dt-dtc2;
		

		  dx_old_j=x_old_j-x_eq[j];
		  dx_old_i=x_old_i-x_eq[j-1];
		   
		  v_old_i=v_old_i*cos(w_p*dt_f)-w_p*dx_old_i*sin(w_p*dt_f);	
		  v_old_j=v_old_j*cos(w_p*dt_f)-w_p*dx_old_j*sin(w_p*dt_f);

		  x_old_i=x_old_i+v_old_i*sin(w_p*dt_f)-dx_old_i*(1-cos(w_p*dt_f));
		  x_old_j=x_old_j+v_old_j*sin(w_p*dt_f)-dx_old_j*(1-cos(w_p*dt_f));

		  //Aplica-se as correcções
		  x[j-1][0]=x_old_i;
		  x[j][0]=x_old_j;
		  v[j-1][0]=v_old_i;
		  v[j][0]=v_old_j;
		  
		  
		}
	   
	    }
	  
	}  
    

      for (Int_t j=0;j<N;j++)
	{
	  //Condições de Fronteira periódicas
	  if (x[j][0]<0)
	    {
	      // cout<< x[j][i]<<endl;;
	      x[j][0]=length;
	      x_eq[j]=x_eq[j] + length;
	      trigger_1=1;
	      //   cout << trigger<<endl;
	    }

	  if (x[j][0]>length)
	    {
	      x[j][0]=0;
	      x_eq[j]=x_eq[j]-length;
	      trigger_2=1;
	    }
	}
      //Reordenação dos vectores consoante a sua posição relativa
 
      //Bubblesort
      if (trigger_1==1 ||trigger_2==1)
	{
	  for(Int_t i=0;i<N-1;i++)
	    {
	      for(Int_t j=i+1;j<N;j++)
		{
		
		  if(x[j][0]<x[i][0])
		    {
		      x_old=x[i][0];
		      x[i][0]=x[j][0];
		      x[j][0]=x_old;

		      v_old=v[i][0];
		      v[i][0]=v[j][0];
		      v[j][0]=v_old;

		      x_eq_old=x_eq[i];
		      x_eq[i]=x_eq[j];
		      x_eq[j]=x_eq_old;
		    
		    }
		}
	    }
        
	  trigger_1=0;
	  trigger_2=0;

	 
	}

      for(Int_t i=0;i<N;i++)
	{	 
	  file << x[i][0] << " " << v[i][0] << " ";
	}
      file << endl;
    }
  file.close();





  

  //-----------------------------------------------------------------------------------------------------------------
  //Output
  Double_t **x_out=new Double_t *[N]; 
  Double_t **v_out=new Double_t *[N];
  Double_t *t_out;
  t_out=new Double_t[N_iter];
  ifstream output;
  output.open("output.txt");
  for(int i=0;i<N;i++)
    {
      x_out[i]=new Double_t[N_iter];
      v_out[i]=new Double_t[N_iter];
    }
  Int_t i=0;
  while(!output.eof())
    {
      output >> t_out[i];
      for(int j=0;j<N;j++)
	{
	  output >>x_out[j][i];
	  output >> v_out[j][i];
	}
      i++;
    }
  output.close();
  i=0;




  //---------------------------------------------------------------------------------------------------------------------------------------

  
  //Outro Campo electrico
  Double_t *E;
  Double_t *x_axis;
  E=new Double_t[(N+1)*parts];
  x_axis=new Double_t[(N+1)*parts];
  
  Double_t **x_aux=new Double_t*[N+1];
  for(Int_t i=0;i<N+1;i++)
    {
      x_aux[i]=new Double_t[parts];       
    }
 
  
  Int_t count=0;
  for(Int_t t_in=0;t_in<N_iter;t_in=t_in+N_iter/10)
    {
      Double_t **x_aux=new Double_t*[N+1];
      for(Int_t i=0;i<N+1;i++)
	{
	  x_aux[i]=new Double_t[parts];       
	}

      x_aux[N][parts-1]=length;

      for(Int_t i=0;i<N+1;i++)
	{
	  if (i==0)
	    {
	      for(Int_t j=0;j<parts;j++)
		{
		  x_aux[i][j]=j*x_out[i][t_in]/(parts);
	 
		}
	    }
	  if(i>0 &&i<N)
	    {
	      x_aux[i][0]=x_out[i-1][t_in];
	      for (Int_t j=1;j<parts;j++)
		{
		  x_aux[i][j]=x_out[i-1][t_in]+ j*(x_out[i][t_in]-x_out[i-1][t_in])/parts;
	    
		}

	    }
	  if(i==N)
	    {
	      x_aux[i][0]=x_out[i-1][t_in];
	 
	      for(Int_t j=1;j<parts-1;j++)
		{
		  x_aux[i][j]=x_out[i-1][t_in]+j*(length-x_out[i-1][t_in])/(parts-1);
	  
		}
	    }
	}

      Double_t *E;
      Double_t *x_axis;
      E=new Double_t[(N+1)*parts];
      x_axis=new Double_t[(N+1)*parts];
      E[0]=0;

  
      for(Int_t i=0;i<N+1;i++)
	{
	  if (i>0)
	    E[i*parts]=E[(i-1)*(parts)+parts-1] + e*n_0*(x_aux[i][0]- x_aux[i-1][parts-1])-e;

	  for(Int_t j=0;j<parts-1;j++)
	    {
	

	      E[i*parts+j+1]=E[i*parts+j]+e*n_0*(x_aux[i][j+1]-x_aux[i][j]);
	      x_axis[i*parts+j]=x_aux[i][j];

	     
	    }
	  x_axis[i*parts+parts-1]=x_aux[i][parts-1];
	}

      c3->SetRightMargin(0.09);
      c3->SetLeftMargin(0.15);
      c3->SetBottomMargin(0.15);

 
      c3_1->SetRightMargin(0.09);
      c3_1->SetLeftMargin(0.15);
      c3_1->SetBottomMargin(0.15);
      //Desenho do gráfico
      if (count<4)
      c3->cd(count+1);

      if (count>=4)
	c3_1->cd(count-3);

 
      E_gr[count]=new TGraph((N+1)*parts,x_axis,E);
      E_gr[count]->SetMarkerStyle(1);
  
 
      E_gr[count]->SetTitle(Form("Electric field/position at time t= %g",t_in*dt));
      E_gr[count]->GetXaxis()->SetTitle("Position [L]");
    
      // E_gr[count]->GetXaxis()->SetLabelSize(12.5);
      E_gr[count]->GetYaxis()->SetTitle("Energy [M V^{2}]");
      // E_gr[count]->GetYaxis()->SetLabelSize(12.5);
      E_gr[count]->GetYaxis()->SetTitleOffset(1);
      auto *axis_1 = E_gr[count]->GetXaxis();
      axis_1->SetLimits(0,length);
      E_gr[count]->Draw();
      count++;
    }

 
  
	
  //----------------------------------------------------------------------------------------------------------------------------------------
  //Gráficos
  c1->cd(1);
  for(int i=0;i<N;i++)
    {
      gr[i]=new TGraph(N_iter,x_out[i],t_out);

      auto *axis = gr[i]->GetXaxis();
      gr[i]->SetMarkerStyle(1);
 
      axis->SetLimits(0,length);
      //   axis->SetLimits(4500,6000);  // along X
    
      /*
	gr[i]->GetHistogram()->SetMaximum(10);   // along          
	gr[i]->GetHistogram()->SetMinimum(-10);  //   Y     
      */
      
      if (i==0)
	{
	  gr[i]->Draw("ap");
	  gr[i]->SetTitle("Position of the sheets along time");
	  gr[i]->GetXaxis()->SetTitle("Position [L]");
	  gr[i]->GetYaxis()->SetTitle("Time [t]");
	}
      else  gr[i]->Draw("p1 same");
    }


  //-----------------------------------------------------------------------------------------------------------------


  
  //-----------------------------------------------------------------------------------------------------------------
  //Integração trapézio
  Double_t *Int_V;
  Int_V=new Double_t[N];
  Double_t Int_V_aux;
  for (Int_t i=0;i<N;i++)
    {
      Int_V_aux=0;
      for(Int_t j=0;j<N_iter;j++)
	{
	  Int_V_aux+=(v_out[i][j+1]+v_out[i][j]);
	}
      Int_V_aux=(Int_V_aux/(2*N_iter));
      Int_V[i]=Int_V_aux;
      h1[0]->Fill(Int_V[i]);
      h1[0]->GetYaxis()->SetTitle("Counts");
      // h1[0]->GetXaxis()->SetTitle("Velocity [L t^(-1)]");
       h1[0]->GetXaxis()->SetTitle("Velocity [L t^{-1}]");

      c2->cd(1);
      h1[0]->Draw();
	
    }

  //----------------------------------------------------------------------------------------------------------------
   //Energia
   c4->cd(1);
  Double_t *Energy;
  Double_t Energy_aux;
  Energy=new Double_t[N_iter];
  
  
    for(Int_t i=0;i<N_iter;i++)
      {
	Energy_aux=0;
	for(Int_t j=0;j<N;j++)
	  {
	    Double_t dx=x_out[j][i]-x_eq[j];
	    Energy_aux=Energy_aux+m*(pow(v_out[j][i],2))+(pow(w_p,2))*m*(pow(dx,2));        
	  }
	Energy[i]=Energy_aux/2;
      }

   c4->SetRightMargin(0.09);
      c4->SetLeftMargin(0.15);
      c4->SetBottomMargin(0.15);
      
    
    Energy_gr=new TGraph(N_iter,t_out,Energy);
 auto *axis2 = Energy_gr->GetXaxis();
 axis2->SetLimits(0,N_iter*dt);
   
  Energy_gr->SetTitle("Energy/Time");
  Energy_gr->GetXaxis()->SetTitle("Time [t]");
  Energy_gr->GetYaxis()->SetTitle("Energy [M V^2]");
  Energy_gr->SetMarkerStyle(1);
   

  Energy_gr->Draw();

    //-----------------------------------------------------------------------------------------------------------------
  //Delete
  
   delete Energy;
  delete E;
  delete x_axis;
  delete[] x_out;
  delete[] v_out;
  delete[] x_aux;
  delete t_out;
  delete x_eq;
  delete[] x;
  delete[] v;
  // delete c1;
  // delete c2;
  // delete h1[1];
}
