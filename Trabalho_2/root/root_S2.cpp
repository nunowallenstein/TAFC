void root_S2()
{
 
  static const Double_t length=10;
  static const Int_t N=100;
  static const Int_t N_iter=4000;
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
  h1[0]=new TH1D("h01","Velocities histogram",1000,-3,3);

 
  //Canvas
  TCanvas *c1 = new TCanvas("c1","Positions", 750, 650);
  TCanvas *c2 = new TCanvas("c2","Velocities", 750, 650);

  c1->Divide(1,1);
  c2->Divide(1,1);

  //--------------------------------------------------------------------------------------------------------------
  //Declaração dos parametros

  //Pos Equilibrio
  Double_t *x_eq;
  x_eq=new Double_t[N];

  //Posição,array 2d para a posição e tempo
  Double_t *x;
  x=new Double_t[N];
  //Velocidade array 2d para a posição e tempo
  Double_t *v;
  v=new Double_t[N];
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
 
      x_eq[i]=(length/(2*N))+ (i * length/N);
      
      // x[i][0]=x_eq[i]+0.5*pow((+1),i);
      
      x[i]=x_eq[i];      
      v[i]=(-0.5)*pow(-1,i);

      if (i==0)
	v[i]=-1;
      
      file << x[i]<< " " << v[i] << " ";
      //Caso as posições iniciais não estejam na caixa
      if (x[i]>length || x[i]<0)
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
  Double_t v_old,x_eq_old;
  Double_t x_old1,v_old1;

  
  //Variável trigger controla consoante a folha debaixo esteja em cima ou a decima em baixo
  Int_t trigger=0;

  //Crossings
  Int_t trigger_2;

  

  for(int i=1;i<N_iter;i++)
    {
      Double_t taux=i*dt;
      file << taux << " ";
      
      for(int j=0;j<N;j++)
	{
	  x_old1=x[j];
	  v_old1=v[j];
	  
	  dx_old=x_old1-x_eq[j];
	  dv_old=v_old1;	  
	  //Força=-dx
	  dv=dv_old-(m*pow(w_p,2)*dx_old*dt);	  
	  //	  incremento nas velocidades

	  v[j]=dv;
	 	  

	  //incremento nas posições
	  dx=dx_old+dv*dt;
	 
	  x[j]=x_eq[j]+dx;
		  
	 
	  //	 cout << x[j][i] <<endl;   	  
	
      

	  /*
	    for(int j=0;j<N;j++)
	    {
	    dx=(+1)*(x[j][i-1]-x_eq[j]);
	    v[j][i]=v[j][i-1]*cos(w_p*dt)-w_p*dx*sin(w_p)*dt;	
	    x[j][i]=x[j][i-1]+v[j][i-1]*sin(w_p*dt)-dx*(1-cos(w_p*dt));
	  */
	  
	  if ((trigger==1 && x[j]>length) || (trigger==2 && x[j]<0 ))
	    {
	      trigger=0;	  
	    }
	  //Condições de Fronteira periódicas
	  if (x[j]<0 && trigger==0)
	    {
	      // cout<< x[j][i]<<endl;;
	      x[j]=length;
	      x_eq[j]=x_eq[j] + length;
	      trigger=1;
	    
	    }

	  if (x[j]>length && trigger==0)
	    {
	      x[j]=0;
	      x_eq[j]=x_eq[j]-length;
	      trigger=2;
	    }      
	}
      for(int j=0;j<N;j++)
	{  
	   
	  //Crossing
	  if(j>0)
	    {	  
	      if ((x[j-1]>x[j]))
		{
		  //Troca-se as velocidades
		  v_old=v[j-1];
		  v[j-1]=v[j];
		  v[j]=v_old;
		  // cout << v_old << endl;
		  //Troca-se os centros de equilibrio
		  x_eq_old=x_eq[j-1];
		  x_eq[j-1]=x_eq[j];
		  x_eq[j]=x_eq_old;
	    
	      
		}
	    }

	  //Relativamente ao caso das condições de fronteira periódicas
	  if (x[N-1]>x[0] &&(trigger ==1||trigger==2))
	    {	     
	      //Troca-se as velocidades
	      v_old= v[N-1];
	      v[N-1]=v[0];
	      v[0]=v_old;

	      //Troca-se os centros de equilibrio
	      x_eq_old= x_eq[N-1];
	      x_eq[N-1]=x_eq[0];
	      x_eq[0]=x_eq_old;
	    }
	  
	  /*
	    if (x[N-1] >x[0]&&trigger==4)
	    {
	    //Troca-se as velocidades
	    v[N-1]=v_old;
	    v[N-1]=v[0];
	    v[0]=v_old;

	    //Troca-se os centros de equilibrio
	    x_eq[N-1]=x_eq_old;
	    x_eq[N-1]=x_eq[0];
	    x_eq[0]=x_eq_old;


	    }
	  */
	  file << x[j] << " " << v[j] << " ";
	}
	  
      file << endl;
    }
  file.close();


    //-----------------------------------------------------------------------------------------------------------------
  //Delete
  
  delete x_eq;
  delete x;
  delete v;


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
      
       if (i==0)  gr[i]->Draw("ap");
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
  
}
