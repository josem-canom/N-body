// Problema de los N cuerpos
// josem-canom

#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;
// Definimos la constante G
const double G=6.674e-11;

double* creavect(int n)
{
    // Función que crea un vector en la memoria dinámica
    double *B;
    B=new double[n];
    return B;
}
 
double** creamat(int i,int j)
{
    // Función que crea una matriz en la memoria dinámica
    double **A;
    A=new double*[i];
    for (int k=0;k<i;k++)
        A[k]=new double[j];
    return A;
}

void borravect(double *v)
{
	// Función que borra un vector de la memoria dinámica
    delete[] v;
}

void borramat (double **m,int i)
{
	// Función que borra una matriz de la memoria dinámica
    for (int k=0;k<i;k++)
        delete[] m[k];
    delete[] m;
}

double norma(double r[3])
{
    // Función para calcular el módulo de un vector
    double sum=0;
    for (int i=0;i<3;i++)
        sum+=r[i]*r[i];
    return sqrt(sum);
}

double* Ag(double *rvtemp,double M,int in,int jn,double *ag)
{
    // Función para calcular la aceleración gravitatoria que genera el cuerpo
    // in de masa "M", sobre el cuerpo jn (cuya masa es irrelevante), de acuerdo
    // a sus respectivas posiciones, que pasamos a la función a través del vector
	// rvtemp, que, como veremos en la función integradora contiene la posición
	// de todos los cuerpos (además de su velocidad). La función también recibe
	// un puntero *ag, dónde escribirá las tres componentes del vector aceleración
	
    double r21[3];
    // Calculamos el vector que une el cuerpo 1 al 2
    for (int i=0;i<3;i++)
        r21[i]=rvtemp[6*jn+i]-rvtemp[6*in+i];
        
    // Calculamos las tres componentes de la aceleración
    for (int i=0;i<3;i++)
        ag[i]=-G*M*r21[i]/pow(norma(r21),3);
        
    // Devolvemos el vector
    return ag;
}

double* Atg(double *rvtemp,double *m,int jn,int nc,double *atg,double *ag)
{
	// La función Atg calcula la aceleración total que recibe el cuerpo jn
	// a causa de la interación gravitatoria con todos los demás, esta función
	// ha de recibir por tanto la posición de todos los cuerpos (rvtemp), las
	// masas de todos ellos (el vector m), el número de cuerpos (nc), y dos
	// punteros a vectores, uno el puntero auxiliar que utiliza la función Ag,
	// puesto que Atg ha de llamarla para llevar a cabo su cometido y el otro
	// otro puntero auxiliar que de nuevo apunta a un vector de tres componentes
	// dónde la función Atg escribirá su resultado.
	
    for (int i=0;i<3;i++)
        atg[i]=0;
    // Calculamos y sumamos la aceleración resultante de todas las fuerzas
    // gravitatorias que ejercen las demás masas sobre la jn
    for (int i=0;i<nc;i++)
    {
    	// Nos aseguramos de no tener en cuenta al propio cuerpo
        if (i!=jn)
        {
        	// Llamamos a Ag que realizará el cálculo
            ag=Ag(rvtemp,m[i],i,jn,ag);
            
            // Y lo sumamos al vector de la aceleración total
            for (int k=0;k<3;k++)
                atg[k]+=ag[k];
        }
    }
    
    // Devolvemos el vector
    return atg;
}

void RK4EMN(double **rvt,double *m,int nc,int n,double tf,double *atg,double *ag)
{
	// Función que integrará por medio del método de Runge Kutta de orden 4
	// nuestro sistema de ecuaciones, para ello recibe la matriz rvt (de 
	// dimensiones n x (6*nc+1), ya que contiene respectivamente, las tres
	// coordenadas del vector posición (x,y,z) y las del vector velocidad, para
	// cada uno de los nc cuerpos, además de una última columna dónde
	// se almacenará el tiempo transcurrido, las n filas corresponden a que n
	// es el número de pasos de integración que le pasaremos a la función), el
	// vector de las masas, como ya hemos mencionado el número de pasos de 
	// integración, el tiempo final y los vectores auxiliares atg y ag necesarios
	// para llamar a las funciones a las que hemos designado el cálculo de la
	// aceleración gravitatoria.
	
	// Dentro de la función definimos el vector auxiliar der que almacenará las
	// derivadas de cada una de las 6nc funciónes que el sistema de ecuaciones
	// diferenciales tiene, del mismo modo que definimos la matriz que contendrá
	// todas las k auxiliares necesarias para la resolción del método de Runge
	// Kutta de orden 4, cuatro por función luego 4*6nc k's distintas; el vector
	// rvtemp auxiliar que contendrá la fila de la matriz rvt correspondiente
	// al paso de integración actual, y dt, el step de integración
	double der[6*nc],k[6*nc][4],rvtemp[6*nc],dt=tf/n;
	
	// Definimos tiempo inicial = 0
	rvt[0][6*nc]=0;
	
	// Iniciamos el bucle para todos los pasos de integración indicados
    for (int i=0;i<n;i++)
    {
    	// Calculamos la columna del tiempo transcurrido
    	rvt[i+1][6*nc]=(i+1)*dt;
    	for (int h=0;h<6*nc;h++)
    		rvtemp[h]=rvt[i][h];
    		
    	// Calculamos la primera k de cada cuerpo y de cada vector de posición
    	// y velocidad
    	for (int j=0;j<nc;j++)
    	{
    		atg=Atg(rvtemp,m,j,nc,atg,ag);
    		for (int l=0;l<6;l++)
    		{
    			if (l>2)
    				k[6*j+l][0]=atg[l-3];
    			else
    				k[6*j+l][0]=rvtemp[6*j+l+3];
    		}
    	}
    	
    	// Calculamos la segunda...
    	for (int j=0;j<6*nc;j++)
    		rvtemp[j]=rvt[i][j]+dt/2*k[j][0];
    	for (int j=0;j<nc;j++)
    	{
    		atg=Atg(rvtemp,m,j,nc,atg,ag);
    		for (int l=0;l<6;l++)
    		{
    			if (l>2)
    				k[6*j+l][1]=atg[l-3];
    			else
    				k[6*j+l][1]=rvtemp[6*j+l+3];
    		}
    	}
    	
    	// ...la tercera...
    	for (int j=0;j<6*nc;j++)
    		rvtemp[j]=rvt[i][j]+dt/2*k[j][1];
    	for (int j=0;j<nc;j++)
    	{
    		atg=Atg(rvtemp,m,j,nc,atg,ag);
    		for (int l=0;l<6;l++)
    		{
    			if (l>2)
    				k[6*j+l][2]=atg[l-3];
    			else
    				k[6*j+l][2]=rvtemp[6*j+l+3];
    		}
    	}
    	
    	// y la cuarta
    	for (int j=0;j<6*nc;j++)
    		rvtemp[j]=rvt[i][j]+dt*k[j][2];
    	for (int j=0;j<nc;j++)
    	{
    		atg=Atg(rvtemp,m,j,nc,atg,ag);
    		for (int l=0;l<6;l++)
    		{
    			if (l>2)
    				k[6*j+l][3]=atg[l-3];
    			else
    				k[6*j+l][3]=rvtemp[6*j+l+3];
    		}
    	}
    	
    	// Calculamos por fin la derivada, con la fórmula conocida y calculamos
    	// las nuevas posiciones y velocidades de todos los cuerpos transcurrido
    	// el intervalo de tiempo dt
    	for (int j=0;j<6*nc;j++)
    	{
			der[j]=(k[j][0]+2*k[j][1]+2*k[j][2]+k[j][3])/6;
			rvt[i+1][j]=rvt[i][j]+dt*der[j];
		}
	}
}

int main()
{
	// Función main
	cout<<"Problema de los N cuerpos"<<endl<<"Jose Manuel Cano Molina"<<endl<<endl;
	cout<<"Calculando..."<<endl<<endl;
	
	// Definimos el número de cuerpos y el número de pasos de integración,
	// hay que tener en cuenta que ambas cosas incrementan en gran medida
	// el tiempo de ejecución necesario
    int nc=5,n=100000;
    
    // Creamos los vectores auxiliares para las funciones encargadas de calcular
    // la aceleración gravitatoria
    double *atg,*ag;
    atg=creavect(3); ag=creavect(3);
    
    // Definimos el tiempo a integrar (en segundos, el vector de las masas, y la
	// matriz rvt
    double *m,**rvt,tf=12e8;
    
    m=creavect(nc);
    // Definimos las masas de los distintos cuerpos
    m[0]=1*1.998e30;
	m[1]=7.349e22;
	m[2]=3e-6*1.998e30;
	m[3]=1e-3*1.998e30;
	m[4]=3e-4*1.998e30;
	
    rvt=creamat(n+1,6*nc+1);
    // Definimos las condiciones iniciales de todos los cuerpos, de la forma en
    // la que habiamos comentado, componente x del primer cuerpo, componente y,
    // componente z, vx, vy, vz, componente x del segundo cuerpo, ..., etc.
    // en unidades de m y m/s
	rvt[0][0]=0;rvt[0][1]=0;rvt[0][2]=0;
	rvt[0][3]=0;rvt[0][4]=0;rvt[0][5]=0;
	rvt[0][6]=0;rvt[0][7]=1.496e11+384400e3;rvt[0][8]=0;
	rvt[0][9]=-3e4-1027;rvt[0][10]=0;rvt[0][11]=0;
	rvt[0][12]=0;rvt[0][13]=1.496e11;rvt[0][14]=0;
	rvt[0][15]=-3e4;rvt[0][16]=0;rvt[0][17]=0;
	rvt[0][18]=0;rvt[0][19]=5.2*1.496e11;rvt[0][20]=0;
	rvt[0][21]=-13.7e3;rvt[0][22]=0;rvt[0][23]=0;
	rvt[0][24]=0;rvt[0][25]=9.537*1.496e11;rvt[0][26]=0;
	rvt[0][27]=-9.67e3;rvt[0][28]=0;rvt[0][29]=0;
	
	// Llamamos a la función integradora
	RK4EMN(rvt,m,nc,n,tf,atg,ag);
	
    // Abrimos el flujo
    ofstream fout("rvt.txt");
    // Sacamos al fichero los datos de la matriz rvt
    for (int i=0;i<=n;i++)
    {
        for (int j=0;j<(6*nc+1);j++)
            fout<<rvt[i][j]<<" ";
        fout<<endl;
    }
    // Cerramos el flujo
    fout.close();
    
    // Indicamos al usuario que los datos han sido transferidos al fichero
    cout<<"Los datos estan en el fichero ''rvt.txt''."<<endl<<endl;
    
    // Borramos toda la memoria dinámica utilizada
    borravect(atg);borravect(ag);borravect(m);borramat(rvt,n+1);
    
    system("PAUSE");
    return 0;
}
