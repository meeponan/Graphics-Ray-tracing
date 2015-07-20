#include <stdio.h>
#define SSD_UTIL_SOURCE_CODE
#include "SSD_util.h"
#define MAXLINELENGTH 1000
#define MAXLABELLEN      16

struct ssd_keyword keyword_table[] = 
  {{SCREEN_KEY, "screen", 5},
   {COLOR_KEY, "color", 3},
   {LINE_KEY, "line", 1},
   {VERTEX_KEY, "vertex", 3},
   {POLYLINE_KEY, "polyline", 2},
   {CIRCLE_KEY, "circle", 3},
   {ARC_KEY, "arc", 5},
   {SAVE_KEY, "save", 1},
   {TRIANGLE_KEY, "triangle",0},
   {EYE_KEY, "eye",3},
   {GAZE_KEY, "gaze",3},
   {UP_KEY, "upvector",3},
   {ORTHO_KEY, "ortho",4},
   {PERS_KEY, "perspective",3},
   {FLOOR_KEY, "floor",5},
   {AXIS_KEY, "axis",2},
   {IDENTITY_KEY, "identity",0}, 
   {TRANSLATE_KEY, "translate",3},
   {ROTATE_KEY, "rotate",3},
   {SCALE_KEY, "scale",3},
   {MESH_KEY, "mesh",2},
   /*project#3 add*/
   {LIGHT_KEY, "light",3},
   {AMBIENT_KEY, "ambient",3},
   {SHADING_KEY, "shading",1},
   {DIFFUSE_KEY, "diffuse",3},
   {SPECULAR_KEY, "specular",3},
   {ATTENUATION_KEY, "attenuation",3},
   {SPOTLIGHT_KEY, "spotlight",5},
   {ALPHA_KEY, "alpha",1},
   {SPHERE_KEY, "sphere",0},
   
   {-1,"unknown", 0}};
    

int match_Keyword(char *keyword, int *npara)
{
  int i;
  *npara = 0;
  for (i=0; keyword_table[i].key_id != -1; i++) {
    if (strcmp(keyword_table[i].name, keyword) == 0) {
      *npara = keyword_table[i].npara;
      return keyword_table[i].key_id;
    }
  }
  return -1;
}

int readAndParse(FILE *inFilePtr, char *keyword, char *arg0,
    char *arg1, char *arg2, char *arg3, char *arg4)
{
  static char line[MAXLINELENGTH];
  char *ptr = line;

  /* delete prior values */
  keyword[0] = arg0[0] = arg1[0] = arg2[0] = arg3[0] = arg4[0] = '\0';
  if (feof(inFilePtr)) return(0);  
  /* read the line from the SSD file */
  while (1) {

        //if(fgetc(inFilePtr) != NULL)
        //fseek(inFilePtr,0,SEEK_CUR);

    if (fgets(line, MAXLINELENGTH, inFilePtr) == NULL) {
      /* reached end of file */
      return(0);
    }
    
    if (feof(inFilePtr)==0) {
      /* check for line too long (by looking for a \n) */
      if (strchr(line, '\n') == NULL) {
    /* line too long */
    printf("error: line too long\n");
    exit(1);
      }
    }
    if (line[0] != '#') break;
  }
  /*
   * Parse the line. 
   */
   if(line[0] == 32)
    {
  sscanf(line, 
  "%*[\t\n ]%[^\t\n ]%*[\t\n ]%[^\t\n ]%*[\t\n ]%[^\t\n ]%*[\t\n ]%[^\t\n ]%*[\t\n ]%[^\t\n ]%*[\t\n ]%[^\t\n ]%*[\t\n ]",
    keyword, arg0, arg1, arg2, arg3, arg4);
    }
    else
    {
  sscanf(line, 
  "%[^\t\n ]%*[\t\n ]%[^\t\n ]%*[\t\n ]%[^\t\n ]%*[\t\n ]%[^\t\n ]%*[\t\n ]%[^\t\n ]%*[\t\n ]%[^\t\n ]%*[\t\n ]",
    keyword, arg0, arg1, arg2, arg3, arg4);
    }
  return strlen(line);
}

int Read_SSD_Scene(char *fname, SCENE *ascene,CAMERA *acamera, char *saved_fname)
{
  char keyword[MAXLINELENGTH],  arg0[MAXLINELENGTH],
    arg1[MAXLINELENGTH], arg2[MAXLINELENGTH], arg3[MAXLINELENGTH];
  char arg4[MAXLINELENGTH];
    int ident = 0;
  FILE *fp;
  int ntranslate = 0;
  int nrotate = 0;
  int nscale = 0;
  int instr_ind = 0;
  int brk = 0;
  int result_readAndParse = 0;
  int shading = 0;
  double diffuse[3] = {0,0,0};
  double specular[3] = {0,0,0};
  double attenuation[3]={1,0,0};
  double alpha = 1;


  /* We first set all the default values */
  RGB_COLOR fcolor, vcolor;
  int ind, ii, num_ver, key_id, key_id1, npara,mm;

  ascene->screen_w = 600;
  ascene->screen_h = 400;
  /* The default color is white */
  ascene->bcolor.rgba[0] = 1.0; 
  ascene->bcolor.rgba[1] = 1.0;
  ascene->bcolor.rgba[2] = 1.0;
  ascene->bcolor.rgba[3] = 1.0;
  fcolor.rgba[0] = 0.0; fcolor.rgba[1] = 0.0; 
  fcolor.rgba[2] = 0.0; fcolor.rgba[3] = 1.0;
  ascene->nlines = 0;
  ascene->npolylines = 0;
  ascene->ntriangles = 0;
  ascene->nidentities = 0;
  ascene->projectType = 0;
  ascene->isAxis = 0;
  ascene->nlights =0;
  ascene->nsphere = 0;
  ascene->nmesh = 0;


  fp = fopen(fname,"rb");
  if (fp == NULL) {
    fprintf(stderr,"%s:%d: Can not open SSD file <%s>.\n",
            __FILE__, __LINE__, fname);
    return -1;
  }
  while (readAndParse(fp, keyword, arg0, arg1, arg2, arg3, arg4) > 0) {
    if (keyword[0] == '\0') {
      /* We simply all blank lines */
      continue;
    }
    key_id = match_Keyword(keyword, &npara);
    switch(key_id) {
    case LINE_KEY:
      ascene->nlines++;
      break;
    case POLYLINE_KEY:
      ascene->npolylines++;
      break;
    case TRIANGLE_KEY:
      ascene->ntriangles++;
      break;
    /*add*/
    case TRANSLATE_KEY:
      ntranslate++;
      break;
    case ROTATE_KEY:
      nrotate++;
      break;
    case SCALE_KEY:
      nscale++;
      break;
    case MESH_KEY:
      ascene->nmesh++;
      break;
    case IDENTITY_KEY:
      ascene->nidentities++;
      break;
    /*project#3 add*/
    case LIGHT_KEY:
      ascene->nlights++;
      break;
    /*project#4 add*/
    case SPOTLIGHT_KEY:
      ascene->nlights++;
      break;
    case SPHERE_KEY:
      ascene->nsphere++;
    break;

    }
  }
  printf("There are %d lines, %d polylines, and %d triangles in %s.\n", 
     ascene->nlines, ascene->npolylines, ascene->ntriangles,
     fname);
  /* We rewind the file to the very beginning to read the file again */
  rewind(fp);
  ascene->lines = (LINE *)malloc(sizeof(LINE) * ascene->nlines);
  ascene->polylines = (POLYLINE *)malloc(sizeof(POLYLINE) * 
                     ascene->npolylines);
  ascene->triangles = (TRIANGLE *)malloc(sizeof(TRIANGLE) * 
                     ascene->ntriangles);
  ascene->nlines = 0;
  ascene->npolylines = 0;
  ascene->ntriangles = 0;

  /*add*/
  ascene->lights = (LIGHT *)malloc(sizeof(LIGHT) * ascene->nlights);
  ascene->identities = (IDENTITY *)malloc(sizeof(IDENTITY) * ascene->nidentities);
  ascene->translate = (TRANSLATE *)malloc(sizeof(TRANSLATE) * ntranslate);
  ascene->rotate = (ROTATE *)malloc(sizeof(ROTATE) * nrotate);
  ascene->scale = (SCALE *)malloc(sizeof(SCALE) * nscale);
  ascene->mesh = (MESH *)malloc(sizeof(MESH) * ascene->nmesh);
  ascene->sphere = (SPHERE *)malloc(sizeof(SPHERE) * ascene->nsphere);
  ascene->nlights = 0;
  ascene->nidentities = 0;
  ntranslate = 0;
  nrotate = 0;
  nscale = 0;
  int nmesh = 0;
  int nsphere = 0;

  while (readAndParse(fp, keyword, arg0, arg1, arg2, arg3, arg4) > 0) {
    if (keyword[0] == '\0') {
      /* We simply all blank lines */
      continue;
    }
    key_id = match_Keyword(keyword, &npara);
    switch(key_id) {

    case SCREEN_KEY:
     ascene->screen_w = atoi(arg0);
     ascene->screen_h = atoi(arg1);
     ascene->bcolor.rgba[0] = atoi(arg2)/255.0;
     ascene->bcolor.rgba[1]  = atoi(arg3)/255.0;
     ascene->bcolor.rgba[2]  = atoi(arg4)/255.0;
     ascene->bcolor.rgba[3] = 1.0;
     break;

    case COLOR_KEY:
      /* We read the color */
      fcolor.rgba[0] = atoi(arg0)/255.0;  fcolor.rgba[1] = atoi(arg1)/255.0;
      fcolor.rgba[2] = atoi(arg2)/255.0;
      break;
    case LINE_KEY:
      ind = ascene->nlines;
      ascene->lines[ind].width = atof(arg0);
      /* We set the default colors */
      vcolor = fcolor;
      num_ver = 0;
      while (readAndParse(fp, keyword, arg0, arg1, arg2, arg3, arg4) > 0) {
        key_id1 = match_Keyword(keyword, &npara);
        switch(key_id1) {
        case VERTEX_KEY:
          ascene->lines[ind].vertices[num_ver].xyzw[0] = atof(arg0);
      ascene->lines[ind].vertices[num_ver].xyzw[1] = atof(arg1);
          ascene->lines[ind].vertices[num_ver].xyzw[2] = atof(arg2);
      memcpy(ascene->lines[ind].vertices[num_ver].rgba,
         vcolor.rgba, sizeof(vcolor.rgba));
          
#if defined(DEBUG_FLAG)
          printf("Point %d %d with color %6.4f %6.4f %6.4f\n",
                 (int)ascene->lines[ind].vertices[num_ver].xyzw[0],
         (int)ascene->lines[ind].vertices[num_ver].xyzw[1],
         ascene->lines[ind].vertices[num_ver].rgba[0],
         ascene->lines[ind].vertices[num_ver].rgba[1],
         ascene->lines[ind].vertices[num_ver].rgba[2]);
#endif
          num_ver ++;
          break;
        case COLOR_KEY:
      vcolor.rgba[0] = atoi(arg0)/255.0; vcolor.rgba[1] = atoi(arg1)/255.0;
          vcolor.rgba[2] = atoi(arg2)/255.0;
          break;
        //default:
         /* printf("%s:%d Line (%s %s %s %s %s %s) ignored.\n",
                 __FILE__, __LINE__, keyword, arg0, arg1, arg2,
                 arg3, arg4);*/
        }
        if (num_ver == 2) {
          break;
        }
      }
      ascene->nlines++;
      break;
    case POLYLINE_KEY:
      ind = ascene->npolylines;
      ascene->polylines[ind].nvertices = atoi(arg0);
      ascene->polylines[ind].width = atof(arg1);
      ascene->polylines[ind].vertices = 
    (COLOR_VERTEX *)malloc(sizeof(COLOR_VERTEX) *
                   ascene->polylines[ind].nvertices);
      vcolor = fcolor;
      num_ver = 0;
      while (readAndParse(fp, keyword, arg0, arg1, arg2, arg3, arg4) > 0) {
        key_id1 = match_Keyword(keyword, &npara);
        switch(key_id1) {
        case VERTEX_KEY:
          ascene->polylines[ind].vertices[num_ver].xyzw[0] = atof(arg0);
      ascene->polylines[ind].vertices[num_ver].xyzw[1] = atof(arg1);
          ascene->polylines[ind].vertices[num_ver].xyzw[2] = atof(arg2);
      memcpy(ascene->polylines[ind].vertices[num_ver].rgba,
         vcolor.rgba, sizeof(vcolor.rgba));
          
#if defined(DEBUG_FLAG)
          printf("Point %d %d with color %6.4f %6.4f %6.4f\n",
                 (int)ascene->polylines[ind].vertices[num_ver].xyzw[0],
         (int)ascene->polylines[ind].vertices[num_ver].xyzw[1],
         ascene->polylines[ind].vertices[num_ver].rgba[0],
         ascene->polylines[ind].vertices[num_ver].rgba[1],
         ascene->polylines[ind].vertices[num_ver].rgba[2]);
#endif
          num_ver ++;
          break;
        case COLOR_KEY:
      vcolor.rgba[0] = atoi(arg0)/255.0; vcolor.rgba[1] = atoi(arg1)/255.0;
          vcolor.rgba[2] = atoi(arg2)/255.0;
          break;
       // default:
          /*printf("%s:%d Line (%s %s %s %s %s %s) ignored.\n",
                 __FILE__, __LINE__, keyword, arg0, arg1, arg2,
                 arg3, arg4); */
        }
        if (num_ver >= ascene->polylines[ind].nvertices) {
          break;
        }
      }
      ascene->npolylines++;
      break;

    case TRIANGLE_KEY:
      ind = ascene->ntriangles;
      vcolor = fcolor;
      num_ver = 0;
      while (readAndParse(fp, keyword, arg0, arg1, arg2, arg3, arg4) > 0) {
        key_id1 = match_Keyword(keyword, &npara);
        switch(key_id1) {
        case VERTEX_KEY:
          ascene->triangles[ind].vertices[num_ver].xyzw[0] = atof(arg0);
      ascene->triangles[ind].vertices[num_ver].xyzw[1] = atof(arg1);
          ascene->triangles[ind].vertices[num_ver].xyzw[2] = atof(arg2);
      memcpy(ascene->triangles[ind].vertices[num_ver].rgba,
         vcolor.rgba, sizeof(vcolor.rgba));  
#if defined(DEBUG_FLAG)
          printf("Point %d %d with color %6.4f %6.4f %6.4f\n",
                 (int)ascene->triangles[ind].vertices[num_ver].xyzw[0],
         (int)ascene->triangles[ind].vertices[num_ver].xyzw[1],
         ascene->triangles[ind].vertices[num_ver].rgba[0],
         ascene->triangles[ind].vertices[num_ver].rgba[1],
         ascene->triangles[ind].vertices[num_ver].rgba[2]);
#endif
          num_ver ++;
          break;
        case COLOR_KEY:
      vcolor.rgba[0] = atoi(arg0)/255.0; vcolor.rgba[1] = atoi(arg1)/255.0;
          vcolor.rgba[2] = atoi(arg2)/255.0;
          break;
        //default:
         /* printf("%s:%d Line (%s %s %s %s %s %s) ignored.\n",
                 __FILE__, __LINE__, keyword, arg0, arg1, arg2,
                 arg3, arg4);  */
        }
        if (num_ver >= 3) {
          break;
        }
      }
      ascene->ntriangles++;
      break;
    
    /*project#3 add*/
        
    case LIGHT_KEY:
    	{
     ascene->lights[ascene->nlights].light[0] = atof(arg0);
     ascene->lights[ascene->nlights].light[1] = atof(arg1);
     ascene->lights[ascene->nlights].light[2] = atof(arg2);

         int result_readAndParse = readAndParse(fp, keyword, arg0, arg1, arg2, arg3, arg4);
         key_id = match_Keyword(keyword, &npara);
            ascene->lights[ascene->nlights].ndirections = 0;
            while(key_id == VERTEX_KEY)
    		{
                ascene->lights[ascene->nlights].directions[ascene->lights[ascene->nlights].ndirections][0] = atof(arg0);
                ascene->lights[ascene->nlights].directions[ascene->lights[ascene->nlights].ndirections][1] = atof(arg1);
                ascene->lights[ascene->nlights].directions[ascene->lights[ascene->nlights].ndirections++][2] = atof(arg2);
                result_readAndParse = readAndParse(fp, keyword, arg0, arg1, arg2, arg3, arg4);
              key_id = match_Keyword(keyword, &npara);
    		}
            fseek(fp,-result_readAndParse,SEEK_CUR);
            ascene->nlights++;
    	}
     break;

    case AMBIENT_KEY:
     ascene->ambient[0] = atof(arg0);
     ascene->ambient[1] = atof(arg1);
     ascene->ambient[2] = atof(arg2);
     break;

    case SHADING_KEY:
      shading = atof(arg0);
     break;

    case DIFFUSE_KEY:
     diffuse[0] = atof(arg0);
     diffuse[1] = atof(arg1);
     diffuse[2] = atof(arg2);
     break;

    case SPECULAR_KEY:
     specular[0] = atof(arg0);
     specular[1] = atof(arg1);
     specular[2] = atof(arg2);
     specular[3] = atof(arg3);
     break;

    case EYE_KEY:
     acamera->eye.xyzw[0] = atof(arg0);
     acamera->eye.xyzw[1] = atof(arg1);
     acamera->eye.xyzw[2] = atof(arg2);
     break;

    case GAZE_KEY:
     acamera->gaze.xyzw[0] = atof(arg0);
     acamera->gaze.xyzw[1] = atof(arg1);
     acamera->gaze.xyzw[2] = atof(arg2);
     break;

    case UP_KEY:
     acamera->up.xyzw[0] = atof(arg0);
     acamera->up.xyzw[1] = atof(arg1);
     acamera->up.xyzw[2] = atof(arg2);
     break;

    case ORTHO_KEY:
     ascene->ortho.right = atof(arg0);
     ascene->ortho.top = atof(arg1);
     ascene->ortho.near = atof(arg2);
     ascene->ortho.far = atof(arg3);
     break;

    case PERS_KEY:
         ascene->projectType = 1;
     ascene->pers.angle = atof(arg0);
     ascene->pers.near = atof(arg1);
     ascene->pers.far = atof(arg2);
     break;

    case FLOOR_KEY:
     ascene->floor.size = atof(arg0);
     ascene->floor.xmin = atof(arg1);
     ascene->floor.xmax = atof(arg2);
     ascene->floor.ymin = atof(arg3);
     ascene->floor.ymax = atof(arg4);
         ascene->floor.color = fcolor;
     break;

    case AXIS_KEY:
         ascene->isAxis = 1;
     ascene->axis.width = atof(arg0);
     ascene->axis.length = atof(arg1);
     break;
     
     /*project4 add*/
     case ATTENUATION_KEY:
        		{
                 attenuation[0] = atof(arg0);
                 attenuation[1] = atof(arg1);
                 attenuation[2] = atof(arg2);
        		}
            	 break;

                case SPOTLIGHT_KEY:
        		{
                 ascene->lights[ascene->nlights].ltype = 1;
                 ascene->lights[ascene->nlights].light[0] = atof(arg0);
                 ascene->lights[ascene->nlights].light[1] = atof(arg1);
                 ascene->lights[ascene->nlights].light[2] = atof(arg2);
                 ascene->lights[ascene->nlights].cutoff = (atof(arg3)/(double)180) * PI;
                 ascene->lights[ascene->nlights].attenuation[0] = atof(arg4);
                 ascene->lights[ascene->nlights].attenuation[1] = attenuation[0];
                 ascene->lights[ascene->nlights].attenuation[2] = attenuation[1];
                 ascene->lights[ascene->nlights].attenuation[3] = attenuation[2];
         ascene->lights[ascene->nlights].ndirections = 1;

                 int result_readAndParse = readAndParse(fp, keyword, arg0, arg1, arg2, arg3, arg4);
                 key_id = match_Keyword(keyword, &npara);
                 if(key_id == VERTEX_KEY)
        		 {
                        ascene->lights[ascene->nlights].directions[0][0] = atof(arg0);
                        ascene->lights[ascene->nlights].directions[0][1] = atof(arg1);
                        ascene->lights[ascene->nlights].directions[0][2] = atof(arg2);
        		 }
            	 else
                        printf("ERROR in spotlight_key1");

                 result_readAndParse = readAndParse(fp, keyword, arg0, arg1, arg2, arg3, arg4);
                 key_id = match_Keyword(keyword, &npara);
                 if(key_id == VERTEX_KEY)
        		 {

                        ascene->lights[ascene->nlights].spot_location[0] = atof(arg0);
                        ascene->lights[ascene->nlights].spot_location[1] = atof(arg1);
                        ascene->lights[ascene->nlights].spot_location[2] = atof(arg2);
                        printf("WE: %f %f %f\n",ascene->lights[ascene->nlights].spot_location[0],ascene->lights[ascene->nlights].spot_location[1],ascene->lights[ascene->nlights].spot_location[2]);
        		 }
            	 else
                        printf("ERROR in spotlight_key2");
                    ascene->nlights++;
        		}
            	break;
            	
                case ALPHA_KEY:
                alpha = atof(arg0);
            	break;
            	/*project4 end*/



    case IDENTITY_KEY:
            result_readAndParse = readAndParse(fp, keyword, arg0, arg1, arg2, arg3, arg4);
            while (result_readAndParse > 0) {
        key_id1 = match_Keyword(keyword, &npara);

        switch(key_id1) {

        case TRANSLATE_KEY:
                    ascene->translate[ntranslate].xyz[0] = atof(arg0);
                    ascene->translate[ntranslate].xyz[1] = atof(arg1);
                    ascene->translate[ntranslate++].xyz[2] = atof(arg2); 
                    ascene->identities[ascene->nidentities].instr[instr_ind++] = TRANSLATE_KEY;
            		break;

        case ROTATE_KEY:
                    ascene->rotate[nrotate].angle = atof(arg0);
                    ascene->rotate[nrotate].xyz[0] = atof(arg1);
                    ascene->rotate[nrotate].xyz[1] = atof(arg2);
                    ascene->rotate[nrotate++].xyz[2] = atof(arg3); 
                    ascene->identities[ascene->nidentities].instr[instr_ind++] = ROTATE_KEY;
            		break;

        case SCALE_KEY:
                    ascene->scale[nscale].xyz[0] = atof(arg0);
                    ascene->scale[nscale].xyz[1] = atof(arg1);
                    ascene->scale[nscale++].xyz[2] = atof(arg2); 
                    ascene->identities[ascene->nidentities].instr[instr_ind++] = SCALE_KEY;
            		break;
        
        /*PROJECT#3 ADD*/
        case LIGHT_KEY:
        		{
                 ascene->lights[ascene->nlights].light[0] = atof(arg0);
                 ascene->lights[ascene->nlights].light[1] = atof(arg1);
                 ascene->lights[ascene->nlights].light[2] = atof(arg2);

                 result_readAndParse = readAndParse(fp, keyword, arg0, arg1, arg2, arg3, arg4);
                 key_id = match_Keyword(keyword, &npara);
                    ascene->lights[ascene->nlights].ndirections = 0;
                    while(key_id == VERTEX_KEY)
        			{
                        ascene->lights[ascene->nlights].directions[ascene->lights[ascene->nlights].ndirections][0] = atof(arg0);
                        ascene->lights[ascene->nlights].directions[ascene->lights[ascene->nlights].ndirections][1] = atof(arg1);
                        ascene->lights[ascene->nlights].directions[ascene->lights[ascene->nlights].ndirections++][2] = atof(arg2);
                        result_readAndParse = readAndParse(fp, keyword, arg0, arg1, arg2, arg3, arg4);
                        key_id = match_Keyword(keyword, &npara);
        			}
                    fseek(fp,-result_readAndParse,SEEK_CUR);
                    ascene->nlights++;
        		}
            	 break;

         case AMBIENT_KEY:
                 ascene->ambient[0] = atof(arg0);
                 ascene->ambient[1] = atof(arg1);
                 ascene->ambient[2] = atof(arg2);
            	 break;

         case SHADING_KEY:
                  shading = atof(arg0);
            	 break;

         case DIFFUSE_KEY:
                 diffuse[0] = atof(arg0);
                 diffuse[1] = atof(arg1);
                 diffuse[2] = atof(arg2);
            	 break;

         case SPECULAR_KEY:
                 specular[0] = atof(arg0);
                 specular[1] = atof(arg1);
                 specular[2] = atof(arg2);
                 specular[3] = atof(arg3);
            	 break;
            	 
       	 /*project#4 add*/
       	 
       	  case ATTENUATION_KEY:
        		{
                 attenuation[0] = atof(arg0);
                 attenuation[1] = atof(arg1);
                 attenuation[2] = atof(arg2);
        		}
            	 break;

       	 case SPOTLIGHT_KEY:
        		{
                 ascene->lights[ascene->nlights].ltype = 1;
                 ascene->lights[ascene->nlights].light[0] = atof(arg0);
                 ascene->lights[ascene->nlights].light[1] = atof(arg1);
                 ascene->lights[ascene->nlights].light[2] = atof(arg2);
                 ascene->lights[ascene->nlights].cutoff = (atof(arg3)/(double)180) * PI;
                 ascene->lights[ascene->nlights].attenuation[0] = atof(arg4);
                 ascene->lights[ascene->nlights].attenuation[1] = attenuation[0];
                 ascene->lights[ascene->nlights].attenuation[2] = attenuation[1];
                 ascene->lights[ascene->nlights].attenuation[3] = attenuation[2];
                 ascene->lights[ascene->nlights].ndirections = 1;

                 int result_readAndParse = readAndParse(fp, keyword, arg0, arg1, arg2, arg3, arg4);
                 key_id = match_Keyword(keyword, &npara);
                 if(key_id == VERTEX_KEY)
        		 {
                        ascene->lights[ascene->nlights].directions[0][0] = atof(arg0);
                        ascene->lights[ascene->nlights].directions[0][1] = atof(arg1);
                        ascene->lights[ascene->nlights].directions[0][2] = atof(arg2);
        		 }
            	 else
                        printf("ERROR in spotlight_key1");

                 result_readAndParse = readAndParse(fp, keyword, arg0, arg1, arg2, arg3, arg4);
                 key_id = match_Keyword(keyword, &npara);
                 if(key_id == VERTEX_KEY)
        		 {
                        ascene->lights[ascene->nlights].spot_location[0] = atof(arg0);
                        ascene->lights[ascene->nlights].spot_location[1] = atof(arg1);
                        ascene->lights[ascene->nlights].spot_location[2] = atof(arg2);
        		 }
            	 else
                        printf("ERROR in spotlight_key2");

                    ascene->nlights++;
        		}
            	break;
    
                case ALPHA_KEY:
                alpha = atof(arg0);
            	break;
   	     
   	     case SPHERE_KEY:
                    ascene->identities[ascene->nidentities].instr[instr_ind++] = SPHERE_KEY;
                    ascene->sphere[nsphere].diffuse[0] = diffuse[0];
                    ascene->sphere[nsphere].diffuse[1] = diffuse[1];
                    ascene->sphere[nsphere].diffuse[2] = diffuse[2];
                    ascene->sphere[nsphere].specular[0] = specular[0];
                    ascene->sphere[nsphere].specular[1] = specular[1];
                    ascene->sphere[nsphere].specular[2] = specular[2];
                    ascene->sphere[nsphere].specular[3] = specular[3];
                    nsphere++;
            	break;
   	    /*project4 end*/

            	 
         case MESH_KEY:
        			{
                        FILE *off = fopen(arg0,"rb");
                        ascene->mesh[nmesh].width = atof(arg1);
                        ascene->mesh[nmesh].alpha = alpha;
                        ascene->mesh[nmesh].shading = shading;
                        for(mm=0; mm < 3; mm ++){
                        	ascene->mesh[nmesh].diffuse[mm] = diffuse[mm];
                        	ascene->mesh[nmesh].specular[mm] = specular[mm];
                        }
                        ascene->mesh[nmesh].specular[3] = specular[3];
                        readAndParse(off, keyword, arg0, arg1, arg2, arg3, arg4);
                        readAndParse(off, keyword, arg0, arg1, arg2, arg3, arg4);
                        ascene->mesh[nmesh].nvertices = atoi(keyword);
                        ascene->mesh[nmesh].vertices = (COLOR_VERTEX *)malloc(sizeof(COLOR_VERTEX) * ascene->mesh[nmesh].nvertices);
                        ascene->mesh[nmesh].npolygons = atoi(arg0);
                        ascene->mesh[nmesh].polygons = (POLYGON *)malloc(sizeof(POLYGON) * ascene->mesh[nmesh].npolygons);
						for(mm = 0; mm < ascene->mesh[nmesh].nvertices; mm++)
        				{
                            readAndParse(off, keyword, arg0, arg1, arg2, arg3, arg4);
                            ascene->mesh[nmesh].vertices[mm].xyzw[0] = atof(keyword);
                            ascene->mesh[nmesh].vertices[mm].xyzw[1] = atof(arg0);
                            ascene->mesh[nmesh].vertices[mm].xyzw[2] = atof(arg1);
                            ascene->mesh[nmesh].vertices[mm].xyzw[3] = 1;
                            if(arg4[0] == '\0')
        					{
                                ascene->mesh[nmesh].vertices[mm].rgba[0] = fcolor.rgba[0];
                                ascene->mesh[nmesh].vertices[mm].rgba[1] = fcolor.rgba[1]; 
                                ascene->mesh[nmesh].vertices[mm].rgba[2] = fcolor.rgba[2];  						
        					}
            				else
        					{
                                ascene->mesh[nmesh].vertices[mm].rgba[0] = atof(arg2);
                                ascene->mesh[nmesh].vertices[mm].rgba[1] = atof(arg3);
                                ascene->mesh[nmesh].vertices[mm].rgba[2] = atof(arg4);
        					}
        				}       				
                        for(mm = 0; mm < ascene->mesh[nmesh].npolygons; mm++)
        				{
                            readAndParse(off, keyword, arg0, arg1, arg2, arg3, arg4);
                            int index[5] = {atoi(arg0),atoi(arg1),atoi(arg2),atoi(arg3),atoi(arg4)};
                            ascene->mesh[nmesh].polygons[mm].nvertices = atoi(keyword);
                			int nn;
                            for(nn = 0; nn < ascene->mesh[nmesh].polygons[mm].nvertices;nn++)
        					{
                                ascene->mesh[nmesh].polygons[mm].num[nn] = index[nn];
        					}
        				}

                	nmesh++;
                    ascene->identities[ascene->nidentities].instr[instr_ind++] = MESH_KEY;

        			}
            		break;        
					        
        case COLOR_KEY:
       /* We read the color */
                    fcolor.rgba[0] = atoi(arg0)/255.0;  fcolor.rgba[1] = atoi(arg1)/255.0;
                    fcolor.rgba[2] = atoi(arg2)/255.0;
                    break;

        case IDENTITY_KEY:
                case SAVE_KEY:
                case LINE_KEY:
                	brk = 1;
                    fseek(fp,-result_readAndParse,SEEK_CUR);
            		break;
        //        default:
         /* printf("%s:%d Line (%s %s %s %s %s %s) ignored.\n",
                 __FILE__, __LINE__, keyword, arg0, arg1, arg2,
                 arg3, arg4);*/
        		}
                if(brk == 1)
            		break;
                result_readAndParse = readAndParse(fp, keyword, arg0, arg1, arg2, arg3, arg4);
    		}
        ascene->identities[ascene->nidentities].instr_num = instr_ind;
        ascene->nidentities++;
        instr_ind = 0;
        brk = 0;
     break;
   /**/

    case SAVE_KEY:
      strcpy(saved_fname, arg0);
      break;
    //default:
     /* printf("%s:%d Keyword (%s) and the line (%s %s %s %s %s) ignored.\n",
             __FILE__, __LINE__, keyword, arg0, arg1, arg2, arg3, arg4);  */
      
    }
  }
  fclose(fp);
  return 0;
}
