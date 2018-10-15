void root_S2()
{
 
  static const Double_t length=10;
  static const Int_t N=20;
  static const Int_t N_iter=10000;
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
  //------------------------------------------------------------------------------------------------------------
  
  //Gráficos
  TGraph *gr[N];


  //Histograma
  TH1D *h1[1];
  h1[0]=new TH1D("h_1","Velocities histogram",1000,-3,3);

 
  //Canvas
  TCanvas *c1 = new TCanvas("c_pos","Positions", 750, 650);
  TCanvas *c2 = new TCanvas("c_Vel","Velocities", 750, 650);

  c1->Divide(1,1);
  c2->Divide(1,1);

  //--------------------------------------------------------------------------------------------------------------
  //Declaração dos parametros

  //Pos Equilibrio
  Double_t *x_eq;
  x_eq=new Double_t[N];

  //Posição,array 2d para a posição e tempo
  Double_t **x=new Double_t *[N];
  
  //Velocidade array 2d para a posição e tempo
  Double_t **v=new Double_t *[N];
  //Tempo
  //  Double_t *t=new Double_t[N_iter];

 
  

  //-----------------------------------------------------------------------------------------------------------------
  
  // TF1 *fa1 = new TF1("fa1","sin(x)/x",0,10);
  // fa1->Draw();
  
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
      // v[i][0]=2.5*pow(-1,i);
      v[i][0]=0;
       if(i==0)
	{
	    v[i][0]=-5;
 
	}

      // if (i==(N/2-1))
      //	v[i]=2;   
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
	  x[j][1]=x[j][0];
	  v[j][1]=v[j][0];
	 
	  dx_old=(+1)*(x[j][1]-x_eq[j]);
	  if(j>0)
	    {
	      //acesso ao x[j-1] da iteração anterior
	      // x_old_i=(x[j-1]-v[j-1]/tan(w_p*dt)+x_eq[j-1]*((w_p)*(tan(w_p*dt))*sin(w_p*dt)+1-cos(w_p*dt)))/((w_p)*(tan(w_p*dt))*sin(w_p*dt)+(cos(w_p*dt)));
	      x_old_i=x[j-1][1];
	      dx_old_i=x_old_i-x_eq[j-1];									    
	      x_old_j=x[j][1];
	      dx_old_j=x_old_j-x_eq[j];

	      v_old_i=v[j-1][1];
	      // v_old_i=v[j-1]/cos(w_p*dt)+w_p*(dx_old_i)*tan(w_p*dt);
	 
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
  //Delete
  
  delete x_eq;
  delete[] x;
  delete[] v;


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
      
       if (i==0) gr[i]->Draw("ap");
        else  gr[i]->Draw("p1 same");
    }
  //-----------------------------------------------------------------------------------------------------------------
  //Integração trapézio
  Double_t *Int_V;
  Int_V=new Double_t[N];
  for (Int_t i=0;i<N;i++)
    {
      for(Int_t j=0;j<N_iter;j++)
	{
	  Int_V[i]+=v_out[i][j]*dt+(v_out[i][j+1]-v_out[i][j])*dt/2;
	}
      Int_V[i]=Int_V[i]/(N_iter*dt);
      h1[0]->Fill(Int_V[i]);
    }

  c2->cd(1);
  h1[0]->Draw();




  delete[] x_out;
  delete[] v_out;
  delete t_out;
  // delete c1;
  // delete c2;
  // delete h1[1];
}
