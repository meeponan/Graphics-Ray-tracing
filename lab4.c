#include <GL/glut.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "SSD_util.h"
#define MAXLINELENGTH 1000
#define MAXFILELEN    256

typedef struct {
    double p0[4];
    double p1[4];
} Line;

typedef struct{
    double rgba[4];
    double z;
    int * bf;
    int * af;
} ALINK;
typedef struct {
    double rgba[4];
    int * start;
    double z;
} HIDDEN;
typedef struct {
    double t;
    double diffuse[3];
    double specular[4];
    double p[3];
    double normal[3];
} record;

char saved_fname[256];
SCENE thescene;
CAMERA vcamera;
HIDDEN *buffer; 
int initial=0;

void display(void);
int Render_SSD(SCENE *ascene, CAMERA *acamera);

void IdentityMatrix(double matrix[][4]);
void MultiMatrix(double matrix1[][4], double matrix2[][4],double Mresult[][4]);
void CoorMatrix(double matrix[][4], double coordinates[]);
void Matrix_Mcam(double Mcam[4][4],double Mcc[][4],double gaze[],double up[],double eye[]);
void Matrix_Mpers(double Mortho[][4],double Minverse[][4], double a,double n,double f,double w,double h);
void Matrix_Mortho(double Mortho[][4],double Minverse[][4],double r,double t,double n,double f);
void Matrix_Mvp(double Mvp[][4], double w, double h);
void MatrixCalculate(double u[],double v[], double t[],double axis[]);
void Matrix_Mrotate(double Mrotate[][4],double axis[]);
void Matrix_Mangle(double Mangle[][4], double axis[],double angle);
void Matrix_Mrotateback(double Mrotateback[][4],double axis[]);
void CrossVector(double vector_left[], double vector_right[], double vector_result[]);

void Matrix_Minverse(double Minverse[][4], double Mcc[][4], double eye[],double w, double h);
void ShadingRender(double Minverse[][4],double Illuminationcolor[], double v0[],double v1[],double v2[],float c0[],float c1[],float c2[],double d[],double s[],int shading);
void Illumination(double Illuminationcolor[],double normal[],double diffuse[],double specular[],double point[]);
void UnitVector(double vector[],double unitvector[]);
/*project#4 add*/ 
void Illumination2(double Illuminationcolor[], double normal[],double diffuse[],double specular[],double illu_point[],double light[],double direction[]);
int raycolor(double Illuminationcolor[], double e[], double d[], double t0, double t1);
int hit(double e[], double d[], double t0, double t1, record* rec);


void Illumination2(double Illuminationcolor[], double normal[],double diffuse[],double specular[],double illu_point[],double light[],double direction[])
{

    int l,d;
    Illuminationcolor[0] = Illuminationcolor[1] = Illuminationcolor[2] = 0;		
   	/*l*/
    UnitVector(direction,direction);
	
	/*n*l*/
	double n_l;
	n_l = normal[0]*direction[0]+normal[1]*direction[1]+normal[2]*direction[2];
	if(n_l < 0){
	    n_l = 0;
	}
	double eye[3] = {vcamera.eye.xyzw[0]-illu_point[0],vcamera.eye.xyzw[1]-illu_point[1],vcamera.eye.xyzw[2]-illu_point[2]};
	UnitVector(eye,eye);
	double half[3] = {eye[0]+direction[0],eye[1]+direction[1],eye[2]+direction[2]};
	UnitVector(half,half);
	/*n*h*/
	double n_h;
	n_h = normal[0]*half[0]+normal[1]*half[1]+normal[2]*half[2];
	if(n_h < 0){
		n_h = 0; }   		
	/*diffuse*/
	Illuminationcolor[0] += diffuse[0] * light[0] * n_l;
	Illuminationcolor[1] += diffuse[1] * light[1] * n_l;
	Illuminationcolor[2] += diffuse[2] * light[2] * n_l;
	/*specular*/
	Illuminationcolor[0] += specular[0] * light[0] * pow(n_h,specular[3]);
	Illuminationcolor[1] += specular[1] * light[1] * pow(n_h,specular[3]);
	Illuminationcolor[2] += specular[2] * light[2] * pow(n_h,specular[3]);
}

int raycolor(double Illuminationcolor[], double e[], double d[], double t0, double t1)
{
    int bkcolor = 0;
    record rec,srec;
    if(hit(e,d,t0,t1,&rec))
    {
        double color[3];
        color[0] = rec.diffuse[0] * thescene.ambient[0];
        color[1] = rec.diffuse[1] * thescene.ambient[1];
        color[2] = rec.diffuse[2] * thescene.ambient[2];
        int i;
        for(i = 0; i < thescene.nlights; i++)
    	{
            double direction[3] = {thescene.lights[i].directions[0][0],
                                         thescene.lights[i].directions[0][1],
                                         thescene.lights[i].directions[0][2]};
            double light[3] = {thescene.lights[i].light[0],thescene.lights[i].light[1],thescene.lights[i].light[2]};

            if(!hit(rec.p,direction,0.01,t1,&srec))
    		{
                    Illumination2(Illuminationcolor,rec.normal,rec.diffuse,rec.specular,rec.p,light,direction);
                    color[0] += Illuminationcolor[0];
                    color[1] += Illuminationcolor[1];
                    color[2] += Illuminationcolor[2];
    		}
    	}

        Illuminationcolor[0] = color[0];
        Illuminationcolor[1] = color[1];
        Illuminationcolor[2] = color[2];

        if(Illuminationcolor[0]>1)
            Illuminationcolor[0] = 1;
        if(Illuminationcolor[1]>1)
            Illuminationcolor[1] = 1;
        if(Illuminationcolor[2]>1)
            Illuminationcolor[2] = 1;
    }
    else
    {
        bkcolor = 1;
    }
    return bkcolor;
}

int hit(double e[], double d[], double t0, double t1, record* rec)
{
    int hit = 0;
    int uu,vv,k,m,n,j;
    rec->t  = t1;
    double t,tt1,tt2;
    
    /*intersect with sphere*/
    for(k = 0; k < thescene.nsphere; k++)
    {
        double es[4] = {e[0],e[1],e[2],e[3]};
        CoorMatrix(thescene.sphere[k].matrix,es);

        double ds[4] = {d[0],d[1],d[2],d[3]};
        CoorMatrix(thescene.sphere[k].matrix,ds);

        double A = ds[0]*ds[0] + ds[1]*ds[1] + ds[2]*ds[2];
        double B = 2*ds[0]*es[0] + 2*ds[1]*es[1] + 2*ds[2]*es[2];
        double C = es[0]*es[0] + es[1]*es[1] + es[2]*es[2] - 1;

        if(B*B - 4*A*C >=0)
    	{
            /*computer t*/
            double t;
            if(B*B - 4*A*C ==0)
    		{
                t = -B/(2*A);
    		}
        	else
    		{
                double t1,t2;
                t1 = (-B + sqrt(B*B - 4*A*C))/(2*A);
                t2 = (-B - sqrt(B*B - 4*A*C))/(2*A);
                if(t1 <= t2)
            		t = t1;
        		else
            		t = t2;
    		}

            if(t < rec->t && t >= t0)
    		{
                double normal[4];
                normal[0] = es[0] + t * ds[0];
                normal[1] = es[1] + t * ds[1];
                normal[2] = es[2] + t * ds[2];
                normal[3] = 1;
                double Mtranspose[4][4];
                for(uu = 0; uu <= 3; uu++ ){
                	for (vv = 0; vv <= 3; vv++){
	                	Mtranspose[uu][vv]=thescene.sphere[k].matrix[vv][uu];
	                }
                }
                CoorMatrix(Mtranspose,normal);
                UnitVector(normal,normal);
                double p[3];
                p[0] = e[0] + t * d[0];
                p[1] = e[1] + t * d[1];
                p[2] = e[2] + t * d[2];
            	hit = 1;
                rec->t = t;
                memcpy(rec->diffuse,thescene.sphere[k].diffuse,sizeof(double)*3);
                memcpy(rec->specular,thescene.sphere[k].specular,sizeof(double)*4);
                memcpy(rec->p,p,sizeof(double)*3);
                memcpy(rec->normal,normal,sizeof(double)*3);

    		}
    	}
    }

    /*intersect with mesh*/
    for(m = 0; m < thescene.nmesh; m++)
    {
        for(n = 0; n < thescene.mesh[m].npolygons; n++)
    	{
            double aa[3] = {thescene.mesh[m].vertices[thescene.mesh[m].polygons[n].num[0]].xyzw[0] - thescene.mesh[m].vertices[thescene.mesh[m].polygons[n].num[1]].xyzw[0],
                       thescene.mesh[m].vertices[thescene.mesh[m].polygons[n].num[0]].xyzw[1] - thescene.mesh[m].vertices[thescene.mesh[m].polygons[n].num[1]].xyzw[1],
                       thescene.mesh[m].vertices[thescene.mesh[m].polygons[n].num[0]].xyzw[2] - thescene.mesh[m].vertices[thescene.mesh[m].polygons[n].num[1]].xyzw[2]};

            double bb[3] = {thescene.mesh[m].vertices[thescene.mesh[m].polygons[n].num[0]].xyzw[0] - thescene.mesh[m].vertices[thescene.mesh[m].polygons[n].num[2]].xyzw[0],
                       thescene.mesh[m].vertices[thescene.mesh[m].polygons[n].num[0]].xyzw[1] - thescene.mesh[m].vertices[thescene.mesh[m].polygons[n].num[2]].xyzw[1],
                       thescene.mesh[m].vertices[thescene.mesh[m].polygons[n].num[0]].xyzw[2] - thescene.mesh[m].vertices[thescene.mesh[m].polygons[n].num[2]].xyzw[2]};
            
            double cc[3] = {thescene.mesh[m].vertices[thescene.mesh[m].polygons[n].num[0]].xyzw[0] - e[0],
                       thescene.mesh[m].vertices[thescene.mesh[m].polygons[n].num[0]].xyzw[1] - e[1],
                       thescene.mesh[m].vertices[thescene.mesh[m].polygons[n].num[0]].xyzw[2] - e[2]};
            	
            double M = aa[0]*(bb[1]*d[2]-d[1]*bb[2]) + aa[1]*(d[0]*bb[2]-bb[0]*d[2]) + aa[2]*(bb[0]*d[1]-bb[1]*d[0]);

            if(M != 0)
    		{
                double t = -(bb[2]*(aa[0]*cc[1]-cc[0]*aa[1]) + bb[1]*(cc[0]*aa[2]-aa[0]*cc[2]) + bb[0]*(aa[1]*cc[2]-cc[1]*aa[2]))/M;
                if(t >= t0 && t < rec->t)
        		{
                    double gamma = (d[2]*(aa[0]*cc[1]-cc[0]*aa[1]) + d[1]*(cc[0]*aa[2]-aa[0]*cc[2]) + d[0]*(aa[1]*cc[2]-cc[1]*aa[2]))/M;
                    if(gamma >= 0 && gamma <= 1)
        			{
                        double beta = (cc[0]*(bb[1]*d[2]-d[1]*bb[2]) + cc[1]*(d[0]*bb[2]-bb[0]*d[2]) + cc[2]*(bb[0]*d[1]-bb[1]*d[0]))/M;
                        if(beta >= 0 && beta <= (1-gamma))
        				{
                    		double p[3];
                            p[0] = e[0] + t * d[0];
                            p[1] = e[1] + t * d[1];
                            p[2] = e[2] + t * d[2];

                			hit = 1;
                    		rec->t = t;
                            memcpy(rec->diffuse,thescene.mesh[m].diffuse,sizeof(double)*3);
                            memcpy(rec->specular,thescene.mesh[m].specular,sizeof(double)*4);
                            memcpy(rec->p,p,sizeof(double)*3);
                            memcpy(rec->normal,thescene.mesh[m].polygons[n].normal,sizeof(double)*3);
        				}
        			}
        		}
    		}
    	}
    	
    }
    return hit;    
}

int max(double a, double b, double c){
	int result=a;
	if(b>result){
		result=b;
	}
	if(c>result){
		result=c;
	}
	return result;
}

int min(double a, double b, double c){
	int result=a;
	if(b<result){
		result=b;
	}
	if(c<result){
		result=c;
	}
	return result;
}


void ShadingRender(double Minverse[][4],double Illuminationcolor[], double v0[],double v1[],double v2[],float c0[],float c1[],float c2[],double d[],double s[],int shading)
{
	int xmint=min(v0[0],v1[0],v2[0]);
	int xmaxt=max(v0[0],v1[0],v2[0]);
	int ymint=min(v0[1],v1[1],v2[1]);
	int ymaxt=max(v0[1],v1[1],v2[1]);

    double line01=(v0[1]-v1[1])*v2[0] + (v1[0] - v0[0])*v2[1] + v0[0]*v1[1] - v1[0]*v0[1];    
    double line02=(v0[1]-v2[1])*v1[0] + (v2[0] - v0[0])*v1[1] + v0[0]*v2[1] - v2[0]*v0[1];
	double line12=(v1[1]-v2[1])*v0[0] + (v2[0] - v1[0])*v0[1] + v1[0]*v2[1] - v2[0]*v1[1];

    double x_incr_a = (v1[1]-v2[1])/line12;
    double x_incr_b = (v0[1]-v2[1])/line02;
    double x_incr_r = (v0[1]-v1[1])/line01;

    double y_incr_a = (v2[0]-v1[0])/line12;
    double y_incr_b = (v2[0]-v0[0])/line02;
    double y_incr_r = (v1[0]-v0[0])/line01;

    double alpha0=((v1[1]-v2[1])*xmint + (v2[0] - v1[0])*ymint + v1[0]*v2[1] - v2[0]*v1[1])/line12;
    double beta0=((v0[1]-v2[1])*xmint + (v2[0] - v0[0])*ymint + v0[0]*v2[1] - v2[0]*v0[1])/line02;
    double gamma0=((v0[1]-v1[1])*xmint + (v1[0] - v0[0])*ymint + v0[0]*v1[1] - v1[0]*v0[1])/line01;

    /* (-1,-1) or (-2,-1) */
    double f12=(v1[1]-v2[1])*(-1) + (v2[0] - v1[0])*(-1) + v1[0]*v2[1] - v2[0]*v1[1];
    double f02=(v0[1]-v2[1])*(-1) + (v2[0] - v0[0])*(-1) + v0[0]*v2[1] - v2[0]*v0[1];
    double f01=(v0[1]-v1[1])*(-1) + (v1[0] - v0[0])*(-1) + v0[0]*v1[1] - v1[0]*v0[1];

    if(f12 == 0)
            f12=(v1[1]-v2[1])*(-2) + (v2[0] - v1[0])*(-1) + v1[0]*v2[1] - v2[0]*v1[1];      
    if(f02 == 0)
            f02=(v0[1]-v2[1])*(-2) + (v2[0] - v0[0])*(-1) + v0[0]*v2[1] - v2[0]*v0[1];
    if(f01 == 0)
            f01=(v0[1]-v1[1])*(-2) + (v1[0] - v0[0])*(-1) + v0[0]*v1[1] - v1[0]*v0[1];

    int x,y;

    double alpha,beta,gamma;	

    for (y = ymint; y <= ymaxt; y++)
    {

        alpha = alpha0;
        beta = beta0;
        gamma = gamma0;

        for (x = xmint; x <= xmaxt; x++)
    	{
            if(alpha >= 0 && beta >= 0 && gamma >= 0)
    		{
                if( (alpha > 0 || line12 * f12 >0) && (beta > 0 || line02 * f02 >0) && (gamma > 0 || line01 * f01 > 0))
              {				
                    double z = alpha * v0[2] + beta * v1[2] + gamma * v2[2];
                    if(z < buffer[y*thescene.screen_w+x].z)
        			{
                        if(shading == 0)
        				{
                            buffer[y*thescene.screen_w+x].rgba[0] = Illuminationcolor[0];
                            buffer[y*thescene.screen_w+x].rgba[1] = Illuminationcolor[1];	
                            buffer[y*thescene.screen_w+x].rgba[2] = Illuminationcolor[2];		
                            buffer[y*thescene.screen_w+x].z = z;
        				}
                        else if(shading == 1)
        				{
                            buffer[y*thescene.screen_w+x].rgba[0] = alpha * c0[0] + beta * c1[0] + gamma * c2[0];
                            buffer[y*thescene.screen_w+x].rgba[1] = alpha * c0[1] + beta * c1[1] + gamma * c2[1];	
                            buffer[y*thescene.screen_w+x].rgba[2] = alpha * c0[2] + beta * c1[2] + gamma * c2[2];		
                            buffer[y*thescene.screen_w+x].z = z;
        				}
            			else
        				{	
                        	double normal[3];
                            normal[0] = alpha * c0[0] + beta * c1[0] + gamma * c2[0];
                            normal[1] = alpha * c0[1] + beta * c1[1] + gamma * c2[1];	
                            normal[2] = alpha * c0[2] + beta * c1[2] + gamma * c2[2];
                            UnitVector(normal,normal);

                            double w = alpha * v0[3] + beta * v1[3] + gamma * v2[3];
                            double p[4] = {x*w,y*w,z,w};
                            CoorMatrix(Minverse,p);
                            Illumination(Illuminationcolor,normal,d,s,p);
                            buffer[y*thescene.screen_w+x].rgba[0] = Illuminationcolor[0];
                            buffer[y*thescene.screen_w+x].rgba[1] = Illuminationcolor[1];	
                            buffer[y*thescene.screen_w+x].rgba[2] = Illuminationcolor[2];		
                            buffer[y*thescene.screen_w+x].z = z;
        				}
        			}
        		}	
    		}

            alpha += x_incr_a;
            beta += x_incr_b;
            gamma += x_incr_r;
    		
    	}

        alpha0 += y_incr_a;
        beta0 += y_incr_b;
        gamma0 += y_incr_r;
    }

}

void Illumination(double Illuminationcolor[],double normal[],double diffuse[],double specular[],double point[])
{

    int mm,nn;
    Illuminationcolor[0] = diffuse[0] * thescene.ambient[0];
    Illuminationcolor[1] = diffuse[1] * thescene.ambient[1];
    Illuminationcolor[2] = diffuse[2] * thescene.ambient[2];

    for(mm = 0; mm < thescene.nlights; mm++)
    {
        for(nn = 0; nn < thescene.lights[mm].ndirections; nn++)
    	{
            double direction[3] = {thescene.lights[mm].directions[nn][0],thescene.lights[mm].directions[nn][1],thescene.lights[mm].directions[nn][2]};
            UnitVector(direction,direction);
            double ndirections;
            ndirections = normal[0]*direction[0]+normal[1]*direction[1]+normal[2]*direction[2];
            if(ndirections < 0){
            	ndirections=0;
            }
            double eye[3] = {vcamera.eye.xyzw[0]-point[0],vcamera.eye.xyzw[1]-point[1],vcamera.eye.xyzw[2]-point[2]};
            UnitVector(eye,eye);
            double half[3] = {eye[0]+direction[0],eye[1]+direction[1],eye[2]+direction[2]};
            UnitVector(half,half);
            double nhalf;
            nhalf = normal[0]*half[0]+normal[1]*half[1]+normal[2]*half[2];
            if(nhalf < 0){
            	nhalf=0;
            }
    		
            /*diffuse*/
            Illuminationcolor[0] += diffuse[0] * thescene.lights[mm].light[0] * ndirections;
            Illuminationcolor[1] += diffuse[1] * thescene.lights[mm].light[1] * ndirections;
            Illuminationcolor[2] += diffuse[2] * thescene.lights[mm].light[2] * ndirections;
            /*specular*/
            Illuminationcolor[0] += specular[0] * thescene.lights[mm].light[0] * pow(nhalf,specular[3]);
            Illuminationcolor[1] += specular[1] * thescene.lights[mm].light[1] * pow(nhalf,specular[3]);
            Illuminationcolor[2] += specular[2] * thescene.lights[mm].light[2] * pow(nhalf,specular[3]);
    	}
    }

    if(Illuminationcolor[0]>1)
        Illuminationcolor[0] = 1;
    if(Illuminationcolor[1]>1)
        Illuminationcolor[1] = 1;
    if(Illuminationcolor[2]>1)
        Illuminationcolor[2] = 1;
}


/*make a vector unit*/
void UnitVector(double vector[],double vector_result[])
{
    double vector_length = sqrt(vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2]);    
    vector_result[0] = vector[0]/vector_length;
    vector_result[1] = vector[1]/vector_length;
    vector_result[2] = vector[2]/vector_length;
}

void CrossVector(double vector_left[], double vector_right[], double vector_result[])
{
    vector_result[0] = vector_left[1] * vector_right[2] - vector_left[2] * vector_right[1];
    vector_result[1] = vector_left[2] * vector_right[0] - vector_left[0] * vector_right[2];
    vector_result[2] = vector_left[0] * vector_right[1] - vector_left[1] * vector_right[0];
}


void CoorMatrix(double matrix[][4], double coordinates[])
{
    int i;
    double sum = 0;
    double coor[4];
    for(i = 0; i < 4; i++)
    {
        sum += matrix[i][0] * coordinates[0];
        sum += matrix[i][1] * coordinates[1];
        sum += matrix[i][2] * coordinates[2];
        sum += matrix[i][3] * coordinates[3];
        coor[i] = sum;
        sum = 0;
    }
    for(i = 0; i < 4; i++)
    {
        coordinates[i] = coor[i];
    }
}

void IdentityMatrix(double matrix[][4])
{
    int i,j;
    for(i = 0; i < 4; i++)
    {
        for(j = 0; j < 4; j++)
        {
            if(i == j)
                matrix[i][j] = 1;
            else
                matrix[i][j] = 0;
        }
    }
}

void MultiMatrix(double matrix1[][4], double matrix2[][4],double Mresult[][4])
{
    double sum = 0;
    double mresult[4][4];
    int i,j;

    for(i = 0; i < 4; i++)
    {
        for(j = 0; j < 4; j++)
        {
            sum += matrix1[i][0] * matrix2[0][j];
            sum += matrix1[i][1] * matrix2[1][j];
            sum += matrix1[i][2] * matrix2[2][j];
            sum += matrix1[i][3] * matrix2[3][j];
            mresult[i][j] = sum;
            sum = 0;
        }
    }
    for(i = 0; i < 4; i++)
    {
        for(j = 0; j < 4; j++)
        {
            Mresult[i][j] = mresult[i][j];
        }
    }
}

void Matrix_Mcam(double Mcam[][4],double Mcc[][4], double gaze[], double up[], double eye[])
{
    double w[3],v[3],u[3];
    UnitVector(gaze,w);
    w[0] = -w[0];
    w[1] = -w[1];
    w[2] = -w[2];
    u[0] = up[1] * w[2] - up[2] * w[1];
    u[1] = up[2] * w[0] - up[0] * w[2];
    u[2] = up[0] * w[1] - up[1] * w[0];
    UnitVector(u,u);
    v[0] = w[1] * u[2] - w[2] * u[1];
    v[1] = w[2] * u[0] - w[0] * u[2];
    v[2] = w[0] * u[1] - w[1] * u[0];
    IdentityMatrix(Mcc);
    IdentityMatrix(Mcam);
    Mcc[0][0] = u[0];
    Mcc[0][1] = u[1];
    Mcc[0][2] = u[2];
    Mcc[1][0] = v[0];
    Mcc[1][1] = v[1];
    Mcc[1][2] = v[2];
    Mcc[2][0] = w[0];
    Mcc[2][1] = w[1];
    Mcc[2][2] = w[2];
    Mcam[0][3] = eye[0];
    Mcam[1][3] = eye[1];
    Mcam[2][3] = eye[2];
    MultiMatrix(Mcc,Mcam,Mcam);
}

void Matrix_Mpers(double Mortho[][4],double Minverse[][4],double a,double n,double f,double w,double h){
    double Mpersinverse[4][4];
    double Morthoinverse[4][4];
    IdentityMatrix(Mortho);
    IdentityMatrix(Mpersinverse);
    IdentityMatrix(Morthoinverse);
    double radian = (a/180.0) * 3.1415926;
    double top = tan(radian/2.0) * (-n);
    double right = (w/h) * top;
    Mortho[0][0] = n/right;
    Mortho[1][1] = n/top;
    Mortho[2][2] = (f + n)/(n - f);
    Mortho[2][3] = (2 * n * f)/(f - n);
    Mortho[3][2] = 1;
    Mortho[3][3] = 0;

    Mpersinverse[0][0] = 1.0/n;
    Mpersinverse[1][1] = 1.0/n;
    Mpersinverse[2][2] = 0;
    Mpersinverse[2][3] = 1;
    Mpersinverse[3][2] = -1.0/(n*f);
    Mpersinverse[3][3] = (n+f)/(n*f);
    	
    Morthoinverse[0][0] = right;
    Morthoinverse[1][1] = top;
    Morthoinverse[2][2] = (n-f)/2.0;
    Morthoinverse[2][3] = (n+f)/2.0;

    MultiMatrix(Minverse,Mpersinverse,Minverse);
    MultiMatrix(Minverse,Morthoinverse,Minverse);
}

void Matrix_Mortho(double Mortho[][4],double Minverse[][4], double r,double t,double n,double f)
{
     IdentityMatrix(Mortho);
     Mortho[0][0] = 1.0/r;
     Mortho[1][1] = 1.0/t;
     Mortho[2][2] = 2.0/(f-n);
     Mortho[2][3] = (f+n)/(n-f);
}

void Matrix_Mvp(double Mvp[][4], double w, double h)
{
    IdentityMatrix(Mvp);
    Mvp[0][0] = w/2.0;
    Mvp[0][3] = (w - 1)/2.0;
    Mvp[1][1] = h/2.0;
    Mvp[1][3] = (h-1)/2.0;
}

void Matrix_Minverse(double Minverse[][4], double Mcc[][4], double eye[],double w, double h){
	double Meye[4][4];
    IdentityMatrix(Meye);
    Meye[0][3] = eye[0];
    Meye[1][3] = eye[1];
    Meye[2][3] = eye[2];
    MultiMatrix(Minverse,Meye,Minverse);
    MultiMatrix(Minverse,Mcc,Minverse);    
    double Mvpinverse[4][4];
    IdentityMatrix(Mvpinverse);
    Mvpinverse[0][0] = 2.0/w;
    Mvpinverse[0][3] = (double)(1.0-w)/w;
    Mvpinverse[1][1] = 2.0/h;
    Mvpinverse[1][3] = (double)(1.0-h)/h;
    MultiMatrix(Minverse,Mvpinverse,Minverse);
}

void MatrixCalculate(double u[],double v[], double t[],double axis[])
{
    UnitVector(axis,axis);
    t[0] = axis[0];
    t[1] = axis[1]+1;
    t[2] = axis[2];
    u[0] = t[1] * axis[2] - t[2] * axis[1];
    u[1] = t[2] * axis[0] - t[0] * axis[2];
    u[2] = t[0] * axis[1] - t[1] * axis[0];
    UnitVector(u,u);
    v[0] = axis[1] * u[2] - axis[2] * u[1];
    v[1] = axis[2] * u[0] - axis[0] * u[2];
    v[2] = axis[0] * u[1] - axis[1] * u[0];
    UnitVector(v,v);
}

void Matrix_Mrotate(double Mrotate[][4],double axis[])
{
    double u[3],v[3],t[3];
    MatrixCalculate(u,v,t,axis);
    Mrotate[0][0] = u[0];
    Mrotate[0][1] = u[1];
    Mrotate[0][2] = u[2];
    Mrotate[1][0] = v[0];
    Mrotate[1][1] = v[1];
    Mrotate[1][2] = v[2];
    Mrotate[2][0] = axis[0];
    Mrotate[2][1] = axis[1];
    Mrotate[2][2] = axis[2];
}

void Matrix_Mangle(double Mangle[][4], double axis[],double angle)
{
    double u[3],v[3],t[3];
    MatrixCalculate(u,v,t,axis);
    double radian = (angle/180.0) * 3.1415926;
    Mangle[0][0] = cos(radian);
    Mangle[0][1] = -sin(radian);
    Mangle[1][0] = sin(radian);
    Mangle[1][1] = cos(radian);
}

void Matrix_Mrotateback(double Mrotateback[][4],double axis[])
{
    double u[3],v[3],t[3];
    MatrixCalculate(u,v,t,axis);
    Mrotateback[0][0] = u[0];
    Mrotateback[0][1] = v[0];
    Mrotateback[0][2] = axis[0];
    Mrotateback[1][0] = u[1];
    Mrotateback[1][1] = v[1];
    Mrotateback[1][2] = axis[1];
    Mrotateback[2][0] = u[2];
    Mrotateback[2][1] = v[2];
    Mrotateback[2][2] = axis[2];
}

void display(void)
{
  Render_SSD(&thescene, &vcamera);
}

int Render_SSD(SCENE *ascene, CAMERA *acamera)
{
  int i,j;
  /* We clear all pixels  */
  glClearColor(ascene->bcolor.rgba[0], ascene->bcolor.rgba[1],
           ascene->bcolor.rgba[2], ascene->bcolor.rgba[3]);
  glClear (GL_COLOR_BUFFER_BIT);  

    buffer = (HIDDEN *)malloc(sizeof(HIDDEN) * ascene->screen_w *ascene->screen_h);

    for(i = 0;i < ascene->screen_w*ascene->screen_h;i++)
    {
        buffer[i].rgba[0] = ascene->bcolor.rgba[0];
        buffer[i].rgba[1] = ascene->bcolor.rgba[1];
        buffer[i].rgba[2] = ascene->bcolor.rgba[2];
        buffer[i].z = 9999;
    }


    double Mfinal[4][4],Mtransform[4][4],Mcam[4][4],Mortho[4][4],Mvp[4][4],Mcc[4][4];
    double Mtranslate[4][4],Mrotate[4][4],Mangle[4][4],Mrotateback[4][4],Mscale[4][4],Mmesh[4][4];
    double gaze[3] = {acamera->gaze.xyzw[0],acamera->gaze.xyzw[1],acamera->gaze.xyzw[2]};
    double up[3] = {acamera->up.xyzw[0],acamera->up.xyzw[1],acamera->up.xyzw[2]};
    double eye[3]= {-acamera->eye.xyzw[0],-acamera->eye.xyzw[1],-acamera->eye.xyzw[2]};
    double Minverse[4][4];
    double Illuminationcolor[3];
    int tt = 0, rr = 0 ,ss = 0, mm = 0, sp = 0,ii,jj;

    /*project#4 add*/
    /*transform to sphere and mesh*/
    if(initial == 0)
    {
        for(i = 0; i < ascene->nidentities; i++)
    	{
            IdentityMatrix(Mtransform);
            for(j = 0; j < ascene->identities[i].instr_num; j++)
    		{
                if(ascene->identities[i].instr[j] == TRANSLATE_KEY)
        		{
                        IdentityMatrix(Mtranslate);
                        Mtranslate[0][3] = ascene->translate[tt].xyz[0];
                        Mtranslate[1][3] = ascene->translate[tt].xyz[1];
                        Mtranslate[2][3] = ascene->translate[tt].xyz[2];
                        MultiMatrix(Mtranslate,Mtransform,Mtransform);
                		tt++;
        		}
                else if(ascene->identities[i].instr[j] == ROTATE_KEY)
        		{
                        double axis[3]={ascene->rotate[rr].xyz[0],ascene->rotate[rr].xyz[1],ascene->rotate[rr].xyz[2]};
                        IdentityMatrix(Mrotate);
						Matrix_Mrotate(Mrotate,axis);
						
						IdentityMatrix(Mangle);
						Matrix_Mangle(Mangle,axis,ascene->rotate[rr].angle);
						MultiMatrix(Mangle,Mrotate,Mrotate);
						
						IdentityMatrix(Mrotateback);
						Matrix_Mrotateback(Mrotateback,axis);
						MultiMatrix(Mrotateback,Mrotate,Mrotate);
						MultiMatrix(Mrotate,Mtransform,Mtransform);
                		rr++;
        			
        		}
                else if(ascene->identities[i].instr[j] == SCALE_KEY)
        		{
                        IdentityMatrix(Mscale);
                        Mscale[0][0] = ascene->scale[ss].xyz[0];
                        Mscale[1][1] = ascene->scale[ss].xyz[1];
                        Mscale[2][2] = ascene->scale[ss].xyz[2];
                        MultiMatrix(Mscale,Mtransform,Mtransform);
                		ss++;
        		}
                else if(ascene->identities[i].instr[j] == SPHERE_KEY)
        		{
                        memcpy(ascene->sphere[sp].matrix,Mtransform,sizeof(double)*16);
                        double tmp[4][4]; 
						double dtmp;
						
  tmp[0][0] = ascene->sphere[sp].matrix[1][1] * ascene->sphere[sp].matrix[2][2] - ascene->sphere[sp].matrix[1][2] * ascene->sphere[sp].matrix[2][1];
  tmp[0][1] = ascene->sphere[sp].matrix[0][2] * ascene->sphere[sp].matrix[2][1] - ascene->sphere[sp].matrix[0][1] * ascene->sphere[sp].matrix[2][2];
  tmp[0][2] = ascene->sphere[sp].matrix[0][1] * ascene->sphere[sp].matrix[1][2] - ascene->sphere[sp].matrix[0][2] * ascene->sphere[sp].matrix[1][1];

  tmp[1][0] = ascene->sphere[sp].matrix[1][2] * ascene->sphere[sp].matrix[2][0] - ascene->sphere[sp].matrix[1][0] * ascene->sphere[sp].matrix[2][2];
  tmp[1][1] = ascene->sphere[sp].matrix[0][0] * ascene->sphere[sp].matrix[2][2] - ascene->sphere[sp].matrix[0][2] * ascene->sphere[sp].matrix[2][0];
  tmp[1][2] = ascene->sphere[sp].matrix[0][2] * ascene->sphere[sp].matrix[1][0] - ascene->sphere[sp].matrix[0][0] * ascene->sphere[sp].matrix[1][2];

  tmp[2][0] = ascene->sphere[sp].matrix[1][0] * ascene->sphere[sp].matrix[2][1] - ascene->sphere[sp].matrix[1][1] * ascene->sphere[sp].matrix[2][0];
  tmp[2][1] = ascene->sphere[sp].matrix[0][1] * ascene->sphere[sp].matrix[2][0] - ascene->sphere[sp].matrix[0][0] * ascene->sphere[sp].matrix[2][1];
  tmp[2][2] = ascene->sphere[sp].matrix[0][0] * ascene->sphere[sp].matrix[1][1] - ascene->sphere[sp].matrix[0][1] * ascene->sphere[sp].matrix[1][0];
  
  dtmp = ascene->sphere[sp].matrix[0][0] * tmp[0][0] + ascene->sphere[sp].matrix[1][0] * tmp[0][1] +
    ascene->sphere[sp].matrix[2][0] * tmp[0][2];
                       
					   dtmp = 1.0/dtmp;
					   	for (ii=0; ii < 3; ii++) {
						   	for (jj=0; jj < 3; jj++) {
							   tmp[ii][jj] *= dtmp;
							   }
						   }
					   	for (jj=0; jj < 3; jj++) {
						   tmp[3][jj] = 0.0;
					   }tmp[3][3] = 1.0;
					   	for (ii=0; ii < 3; ii++) {
						   tmp[ii][3] = 0;
						   	for (jj=0; jj < 3; jj++) {
							   tmp[ii][3] -=  tmp[ii][jj] * ascene->sphere[sp].matrix[jj][3];
							   }
						   }
					   memcpy(ascene->sphere[sp].matrix,tmp,sizeof(double)*16);
					   sp++;
        		}
                else if(ascene->identities[i].instr[j] == MESH_KEY)
        		{
            			int k;
                        for(k = 0; k < ascene->mesh[mm].nvertices; k++)
        				{
                            CoorMatrix(Mtransform,ascene->mesh[mm].vertices[k].xyzw);
        				}

                        for(k = 0; k < ascene->mesh[mm].npolygons; k++)
        				{
                            double normal[3];                             
                            double ab[3];
    ab[0] = ascene->mesh[mm].vertices[ascene->mesh[mm].polygons[k].num[1]].xyzw[0] - ascene->mesh[mm].vertices[ascene->mesh[mm].polygons[k].num[0]].xyzw[0];
    ab[1] = ascene->mesh[mm].vertices[ascene->mesh[mm].polygons[k].num[1]].xyzw[1] - ascene->mesh[mm].vertices[ascene->mesh[mm].polygons[k].num[0]].xyzw[1];
    ab[2] = ascene->mesh[mm].vertices[ascene->mesh[mm].polygons[k].num[1]].xyzw[2] - ascene->mesh[mm].vertices[ascene->mesh[mm].polygons[k].num[0]].xyzw[2];
                            double ac[3];
    ac[0] = ascene->mesh[mm].vertices[ascene->mesh[mm].polygons[k].num[2]].xyzw[0] - ascene->mesh[mm].vertices[ascene->mesh[mm].polygons[k].num[0]].xyzw[0];
    ac[1] = ascene->mesh[mm].vertices[ascene->mesh[mm].polygons[k].num[2]].xyzw[1] - ascene->mesh[mm].vertices[ascene->mesh[mm].polygons[k].num[0]].xyzw[1];
    ac[2] = ascene->mesh[mm].vertices[ascene->mesh[mm].polygons[k].num[2]].xyzw[2] - ascene->mesh[mm].vertices[ascene->mesh[mm].polygons[k].num[0]].xyzw[2];
	                        CrossVector(ab,ac,normal);
                            UnitVector(normal,normal); 
							       					
                            memcpy(ascene->mesh[mm].polygons[k].normal,normal,sizeof(double)*3);
        				}
                		mm++;
        		}
        		else
        		{
                        printf("unknown error!!!");
        		}		
    		}
    	}

        initial = 1;
    }


    /*uvw*/
    double w[3],u[3],v[3];

    UnitVector(gaze,w);
    w[0] = -w[0];
    w[1] = -w[1];
    w[2] = -w[2];
    CrossVector(up,w,u);
    UnitVector(u,u);    
    CrossVector(w,u,v); 

    double radian = (ascene->pers.angle/180.0) * 3.1415926;
    double top = tan(radian/2.0) * (-ascene->pers.near);
    double right = ((double)ascene->screen_w/(double)ascene->screen_h) * top;

    /*each pixel in screen do*/
    for(i = 0; i < ascene->screen_w; i++)
    {
        for(j = 0; j < ascene->screen_h; j++)
    	{
            /*screen to eye space*/
            double us = -right + (2*right)*(((double)i+0.5)/ascene->screen_w);
            double vs = -top + (2*top)*(((double)j+0.5)/ascene->screen_h);
            double ws = ascene->pers.near;

            double e[4] = {acamera->eye.xyzw[0],acamera->eye.xyzw[1],acamera->eye.xyzw[2],1};
            double d[4] = {us*u[0]+vs*v[0]+ws*w[0],
                                    	 us*u[1]+vs*v[1]+ws*w[1],
                                         us*u[2]+vs*v[2]+ws*w[2],0};

            if(!raycolor(Illuminationcolor,e,d,0,999))
    		{
                glBegin(GL_POINTS);
                glColor3f(Illuminationcolor[0],Illuminationcolor[1],Illuminationcolor[2]);
                glVertex2i(i,j);
            	glEnd();
    		}

    	}
    }

  glFlush();
  glutSwapBuffers();

  return 0;
}






void init (void)
{
  /* select clearing color  to the specified background  */
  glClearColor(thescene.bcolor.rgba[0], thescene.bcolor.rgba[1], 
           thescene.bcolor.rgba[2], thescene.bcolor.rgba[3]);
  glClear (GL_COLOR_BUFFER_BIT);
  /* initialize viewing values  */
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0.0, thescene.screen_w-1.0, 
      0.0, thescene.screen_h-1.0, -1.0, 1.0);
  
}

int main(int argc, char** argv)
{
  int ii, jj, kk, argc_1,dd;
  char **my_argv;
  char ssd_fname[MAXFILELEN];
  if (argc < 2) {
    printf("%s:%d Usage: %s SSD_file\n", 
       __FILE__, __LINE__, argv[0]);
    return 0;
  }
  strcpy(ssd_fname, argv[1]);
  strcpy(saved_fname,"graphics_tmp.ppm");
  argc_1 = argc - 1;
  my_argv = (char **)malloc(sizeof(char *) * argc);
  my_argv[0] = argv[0];
  for (ii=2; ii <= argc; ii++) {
    my_argv[ii-1] = argv[ii];
  }
  glutInit(&argc_1, my_argv);
  free(my_argv);
  glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB);
  /* Set the default size and background */
  Read_SSD_Scene(ssd_fname, &thescene, &vcamera,saved_fname);
  glutInitWindowSize (thescene.screen_w, thescene.screen_h);
  glutInitWindowPosition (50, 50);
  glutCreateWindow (argv[0]);
  init ();
  glutDisplayFunc(display);  
 // glutMouseFunc(mouse);
 // glutKeyboardFunc(ssd_keyboard);
  glutMainLoop();
  return 0;   /* ANSI C requires main to return int; it will never be 
         reached as glutMainLoop() does not return. */
}

