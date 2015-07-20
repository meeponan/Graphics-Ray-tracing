#if !defined(SSD_UTIL_H_H)
#define SSD_UTIL_H_H
#include <stdlib.h>
#include <string.h>

#define PI 3.1415926
#define SCREEN_KEY   0
#define COLOR_KEY    1
#define LINE_KEY     2
#define VERTEX_KEY   3
#define POLYLINE_KEY 4
#define CIRCLE_KEY   5
#define ARC_KEY      6
#define SAVE_KEY     7
#define TRIANGLE_KEY 8
#define EYE_KEY      9
#define GAZE_KEY    10
#define UP_KEY      11
#define ORTHO_KEY   12
#define PERS_KEY    13
#define FLOOR_KEY   14
#define AXIS_KEY    15
#define IDENTITY_KEY   16
#define TRANSLATE_KEY  17
#define ROTATE_KEY     18
#define SCALE_KEY      19
#define MESH_KEY       20

#define LIGHT_KEY			21	
#define AMBIENT_KEY		22
#define SHADING_KEY		23
#define DIFFUSE_KEY		24
#define SPECULAR_KEY	25
#define ATTENUATION_KEY	26
#define SPOTLIGHT_KEY	27
#define ALPHA_KEY 28
#define SPHERE_KEY 29


struct ssd_keyword {
  /* Keyword table entry to be used for reading SSD */
  int key_id;
  char name[32];
  int npara;
};

typedef struct {
  float rgba[4];
} RGB_COLOR;

typedef struct {
  double xyzw[4];
} VERTEX;

typedef struct {
  double xyzw[4];
  float rgba[4];
} COLOR_VERTEX;

typedef struct {
  double width;
  COLOR_VERTEX vertices[2];
} LINE;

typedef struct {
  double width;
  int    nvertices;
  COLOR_VERTEX *vertices;
} POLYLINE;

typedef struct {
  COLOR_VERTEX vertices[3];
} TRIANGLE;

/*project3 add*/ 
typedef struct {
  double xyz[3];
} TRANSLATE;

typedef struct {
	double angle;
  double xyz[3];
} ROTATE;

typedef struct {
  double xyz[3];
} SCALE;

typedef struct {
  int nvertices;
	int num[50];
	double normal[3];
} POLYGON;

typedef struct {
  double width;
	int shading;
	double diffuse[3];
	double specular[4];
	int nvertices;
	COLOR_VERTEX *vertices;
	int npolygons;
	POLYGON *polygons;
	double alpha;
	
} MESH;

/*project 4 add*/
typedef struct {
double diffuse[3];
double specular[4];
double matrix[4][4];	
} SPHERE;

typedef struct {
  double size;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  RGB_COLOR color;
} FLOOR;

typedef struct {
	double angle;
	double near;
	double far;
} PERS;

typedef struct {
	double right;
	double top;
	double near;
	double far;
} ORTHO;

typedef struct {
	double width;
	double length;
} AXIS;

typedef struct {
	int instr_num;
	int instr[50];
} IDENTITY;

typedef struct {
	int ltype;
	double light[3]; 
	int ndirections; 
	double directions[10][3];	
	double cutoff;
	double attenuation[4];
	double spot_location[3];
} LIGHT;


typedef struct {
  int screen_w, screen_h;
  RGB_COLOR bcolor;
  int projectType;
  ORTHO ortho;
  PERS pers;
  FLOOR floor;
  int isAxis;
  AXIS axis;
  int nidentities;
  IDENTITY *identities;
  TRANSLATE *translate;
  ROTATE *rotate;
  SCALE *scale;
  int nmesh;
  MESH *mesh;
  int nsphere;
  SPHERE * sphere;
  int nlights;
  LIGHT *lights;
  double ambient[3];
  int  nlines; 
  LINE *lines;
  int  npolylines;
  POLYLINE *polylines;
  int  ntriangles; 
  TRIANGLE *triangles;
} SCENE;

typedef struct {
	VERTEX eye;
	VERTEX gaze;
	VERTEX up;
} CAMERA;

#if defined(SSD_UTIL_SOURCE_CODE)
#define EXTERN_FLAG
#else
#define EXTERN_FLAG extern
extern struct ssd_keyword keyword_table[];
#endif

EXTERN_FLAG 
int match_Keyword(char *keyword, int *npara);

EXTERN_FLAG 
int readAndParse(FILE *inFilePtr, char *keyword, char *arg0,
		 char *arg1, char *arg2, char *arg3, char *arg4);
EXTERN_FLAG
int Read_SSD_Scene(char *fname, SCENE *ascene, CAMERA *acamera, char *saved_fname);
#undef  EXTERN_FLAG

#endif
