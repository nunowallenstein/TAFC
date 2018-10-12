
void root_S2()
{
  static const Double_t length=10;
  static const Int_t N=4;
  static const Int_t N_iter=10000;
  static const Double_t dt =0.01;

  //Vamos por numa file os tempos, as posições e as velocidades
  ofstream file;
  file.open("output.txt");
  


  //Aquisição dos argumentos, comprimento da caixa e o número de folhas 

   
  //Posições iniciais
  Double_t x_i=1;
  Double_t v_i=0;

  
    
  /*
    cout << "teste" << endl;
    cout	<< "Nsheet=" << N<< endl;
    cout	<<"length="<<length<<endl ;
  */
   TGraph *gr[N];
  
  Double_t **x_eq=new Double_t * [N];
  Double_t **x=new Double_t* [N];
  Double_t **v=new Double_t* [N];
  Double_t **v_eq=new Double_t* [N];
  Double_t *t=new Double_t[N_iter];

  

  //Inicialização
  file << "0 "; 
  for (int i=0;i<N;i++)
    {
      x_eq[i]=new Double_t[N_iter];
      v_eq[i]=new Double_t[N_iter]; 
      x[i]=new Double_t[N_iter];
      v[i]=new Double_t[N_iter];
      
      x_eq[i][0]=(length/(2*N))+ (i * length/N);
      x[i][0]=i;     
      v[i][0]=v_i;
      v_eq[i][0]=0;
      file << x[i][0] << " " << v[i][0] << " ";
    }
  
  
  file << endl;
  Double_t dx,dv;
  for(int i=1;i<N_iter;i++)
    {
      Double_t taux=i*dt;
      t[i]=taux;
      file << taux << " ";
      
      for(int j=0;j<N;j++)
	{
	  
	  dx=x[j][i-1]-x_eq[j][i-1];
	  // cout << dx<<endl;
	  dv=v[j][i-1]-v_eq[j][i-1];
	  
	  //Força=-dx
	  dv=dv-dx*dt;	  
	   //	  incremento nas velocidades

	  v[j][i]=dv+v_eq[j][i-1];
	  v_eq[j][i]=v_eq[j][i-1];	  

	  //incremento nas posições
	    dx=dx+dv*dt;
	 
	  x[j][i]=x_eq[j][i-1]+dx;
	  x_eq[j][i]=x_eq[j][i-1];	  
	 
	      cout << x[j][i] <<endl;   
	  file << x[j][i] << " " << v[j][i] << " ";	  
	}
      /*
      //Método stormer verlet
      for(j=0;j<N;j++)
	{
 dx_old=x[j][i-1]-x_eq[j][i-1];
 dv_old=v[j][i-1]-v_eq[j][i-1];

 dx=dx_old+dv_old*dt -1/2*dx_old*dt^2;
 dv=dv_old+(1/2)*(-dx_old-dx)*dt;
 x[j][i]=x_eq[j][i]+dx
	}
      */
      file << endl;
    }
  
  //Output
   for(int i=0;i<N;i++)
   {
      gr[i]=new TGraph(N_iter,t,x[i]);
      if (i==0)
      gr[i]->Draw();

      else
	gr[i]->Draw("same");
   }

   
  delete[] x_eq;
  delete[] v_eq;
  delete[] x;
  delete[] v;

  file.close();
  
}
