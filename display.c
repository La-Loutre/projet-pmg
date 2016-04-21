
#define _XOPEN_SOURCE 600

#include "display.h"

#ifdef __APPLE__
#include <OpenGL/gl.h>          // Header File For The OpenGL32 Library
#include <OpenGL/glu.h>         // Header File For The GLu32 Library
#include <GLUT/glut.h>          // Header File For The GLut Library
#else
#define GL_GLEXT_PROTOTYPES
#include <GL/gl.h>          // Header File For The OpenGL32 Library
#include <GL/glu.h>         // Header File For The GLu32 Library
#include <GL/glut.h>          // Header File For The GLut Library
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <sys/time.h>

#define PIX_SURFACE

#ifdef PIX_SURFACE
unsigned flatten_surface = 0;
#endif

GLfloat *s_vbo_vertex;
GLfloat *s_vbo_color;
GLuint *s_vbo_idx;
GLuint s_vi = 0, s_ci = 0;
GLuint s_nbv = 0;
GLuint s_vertices = 0;
GLuint s_indexes = 0;
GLuint s_vbovid, s_vboidx, s_vbocid, s_texid;

unsigned SAND_DIM = 0;
unsigned SAND_MAX_HEIGHT = 0;

struct timeval lastDisplayTimeval ;

#define MAX(a, b) ((a) > (b) ? (a) : (b))

#define winx 800
#define winy 800

#define EXTENT  5.0

float min_ext[3] = {0.0, 0.0, 0.0},
  max_ext[3] = {EXTENT, EXTENT, EXTENT};

GLuint glutWindowHandle = 0;
GLdouble fovy, aspect, near_clip, far_clip;  
                          /* parameters for gluPerspective() */
GLfloat near;

// mouse controls
GLfloat angle = 0.0;
float dis;

float eye[3];             /* position of eye point */
float center[3];          /* position of look reference point */
float up[3];              /* up direction for camera */

unsigned camera_rotate = 0;
unsigned rescale = 0;

int mouse_old_x, mouse_old_y;
int mouse_buttons = 0;
float rotate_x = 0.0, rotate_y = 0.0;
float translate_z = -1.f;

get_func_t get_func_ptr;
compute_func_t compute_func_ptr;

GLfloat *texture = NULL;


void initializeVBO()
{
#ifdef PIX_SURFACE
  s_vbo_vertex = malloc (4* 3 * SAND_DIM * SAND_DIM * sizeof(GLfloat));
  s_vbo_color = malloc (4 * 3 * SAND_DIM * SAND_DIM * sizeof(GLfloat));
  s_vbo_idx = malloc (4 * 4 * (SAND_DIM) * (SAND_DIM) * sizeof(GLuint));
#else
  s_vbo_vertex = malloc (3 * SAND_DIM * SAND_DIM * sizeof(GLfloat));
  s_vbo_color = malloc (3 * SAND_DIM * SAND_DIM * sizeof(GLfloat));
  s_vbo_idx = malloc (4 * (SAND_DIM-1) * (SAND_DIM-1) * sizeof(GLuint));
#endif

  s_vi = 0;
  s_nbv = 0;
  s_vertices = 0;
  s_indexes = 0;
}

void vboFinalize()
{
  free(s_vbo_vertex);
  free(s_vbo_color);
  free(s_vbo_idx);

  glDeleteBuffers(1, &s_vbocid);
  glDeleteBuffers(1, &s_vboidx);
  glDeleteBuffers(1, &s_vbovid);
}

unsigned color_index;

static void add_sand_color (GLfloat x, GLfloat y, GLfloat z)
{
#ifdef PIX_SURFACE
  if (texture) {
    s_vbo_color[s_ci++] = texture[color_index];
    s_vbo_color[s_ci++] = texture[color_index + 1];
    s_vbo_color[s_ci++] = texture[color_index + 2];

  } else {
    float ratio = y / EXTENT;

    if (ratio > 1.0)
      ratio = 1.0;

    s_vbo_color[s_ci++] = ratio;
    s_vbo_color[s_ci++] = 0.0;
    s_vbo_color[s_ci++] = 1.0 - ratio;
  }
#else
  float ratio = y / EXTENT;

  if (ratio > 1.0)
    ratio = 1.0;

  s_vbo_color[s_ci++] = ratio;
  s_vbo_color[s_ci++] = 0.0;
  s_vbo_color[s_ci++] = 1.0 - ratio;
#endif
}

static void add_sand_vertice (GLfloat x, GLfloat y, GLfloat z, int set_color)
{  
  s_vbo_vertex[s_vi++] = x;
#ifdef PIX_SURFACE
  s_vbo_vertex[s_vi++] = (flatten_surface ? 0 : y);
#endif
  s_vbo_vertex[s_vi++] = z;

  s_nbv++;
  if (set_color)
    add_sand_color (x, y, z);
}

static unsigned max_value;

void sand_surface_refresh (int set_color)
{
  s_vi = s_ci = s_nbv = 0;
  color_index = 0;

  // build vertices
  max_value = 0;
  for (int y = 0; y < SAND_DIM; y++)
    for (int x = 0; x < SAND_DIM; x++) {
      GLfloat xv, yv, zv;
      unsigned val = get_func_ptr (y, x) ;

      if (val > max_value)
	max_value = val;

#ifdef PIX_SURFACE
      color_index = 3 * (y * SAND_DIM + x);

      xv = x * EXTENT / (SAND_DIM);
      yv = val * EXTENT / SAND_MAX_HEIGHT;
      zv = y * EXTENT / (SAND_DIM);
      // 1st vertex
      add_sand_vertice (xv, yv, zv, set_color);
      // 2nd vertex
      add_sand_vertice (xv + EXTENT/SAND_DIM, yv, zv, set_color);
      // 3rd vertex
      add_sand_vertice (xv + EXTENT/SAND_DIM, yv, zv + EXTENT/SAND_DIM, set_color);
      // 4th vertex
      add_sand_vertice (xv, yv, zv + EXTENT/SAND_DIM, set_color);
#else
      xv = x * EXTENT / (SAND_DIM-1);
      yv = val * EXTENT / SAND_MAX_HEIGHT;
      zv = y * EXTENT / (SAND_DIM-1);
      add_sand_vertice (xv, yv, zv, set_color);
#endif
    }
}

void sand_surface_build (void)
{
  sand_surface_refresh (1);

  s_vertices = s_nbv;
  printf ("Sable vertices = %d\n", s_vertices);

#ifdef PIX_SURFACE
  for (int y = 0; y < SAND_DIM; y++)
    for (int x = 0; x < SAND_DIM; x++) {
      int base = 4 * (y * SAND_DIM + x);
      // horiz quad
      s_vbo_idx[s_indexes++] = base;
      s_vbo_idx[s_indexes++] = base + 1;
      s_vbo_idx[s_indexes++] = base + 2;
      s_vbo_idx[s_indexes++] = base + 3;
      // 1st vertic quad
      if (x < SAND_DIM - 1) {
	s_vbo_idx[s_indexes++] = base + 1;
	s_vbo_idx[s_indexes++] = base + 4;
	s_vbo_idx[s_indexes++] = base + 7;
	s_vbo_idx[s_indexes++] = base + 2;
      }
      // 2nd vertic quad
      if (y < SAND_DIM - 1) {
	int base2 = base + 4 * SAND_DIM;
	s_vbo_idx[s_indexes++] = base + 3;
	s_vbo_idx[s_indexes++] = base + 2;
	s_vbo_idx[s_indexes++] = base2 + 1;
	s_vbo_idx[s_indexes++] = base2;
      }
    }
#else
  // build quads
  for (int y = 0; y < SAND_DIM - 1; y++)
    for (int x = 0; x < SAND_DIM - 1; x++) {
      s_vbo_idx[s_indexes++] = y * SAND_DIM + x;
      s_vbo_idx[s_indexes++] = y * SAND_DIM + x + 1;
      s_vbo_idx[s_indexes++] = (y + 1) * SAND_DIM + x + 1;
      s_vbo_idx[s_indexes++] = (y + 1) * SAND_DIM + x;
    }
#endif
  printf ("Sable quads edges = %d\n", s_indexes);

}

void buildVBO()
{
  sand_surface_build ();
  
  glGenBuffers(1, &s_vbocid);
  glBindBuffer(GL_ARRAY_BUFFER, s_vbocid);
  glColorPointer(3, GL_FLOAT, 0, 0);
  glBufferData(GL_ARRAY_BUFFER, s_vertices*3*sizeof(float), s_vbo_color, GL_DYNAMIC_DRAW);

  glGenBuffers(1, &s_vbovid);
  glBindBuffer(GL_ARRAY_BUFFER, s_vbovid);
  glVertexPointer(3, GL_FLOAT, 0, 0);
  glBufferData(GL_ARRAY_BUFFER, s_vertices*3*sizeof(float), s_vbo_vertex, GL_DYNAMIC_DRAW);

  glGenBuffers(1, &s_vboidx);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, s_vboidx);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, s_indexes*sizeof(GLuint), s_vbo_idx, GL_STATIC_DRAW);

  glGenTextures(1, &s_texid); /* Texture name generation */
  glBindTexture(GL_TEXTURE_2D, s_texid); /* Binding of texture name */
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
}

static void render_sand_surface (void)
{
  glEnableClientState(GL_COLOR_ARRAY);
  glEnableClientState(GL_VERTEX_ARRAY);
  
  glDrawElements(GL_QUADS, s_indexes, GL_UNSIGNED_INT, 0);
  
  glDisableClientState(GL_COLOR_ARRAY);
  glDisableClientState(GL_VERTEX_ARRAY);
}

static void reshape (int w, int h)
{
  // set the GL viewport to match the full size of the window
  glViewport(0, 0, (GLsizei)w, (GLsizei)h);
  aspect = w/(float)h;
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(fovy,aspect,near_clip,far_clip);
  glMatrixMode(GL_MODELVIEW);
}

static unsigned true_redisplay = 1;

void display_draw_scene() {

  if(camera_rotate) {
    float p, q;
    q = rotate_y + 180.0; q = q/360.0 * 2 * M_PI;
    rotate_y -= 0.25;
    while(rotate_y < 0.0)
      rotate_y += 360.0;

    p = rotate_y + 180.0; p = p/360.0 * 2 * M_PI;
    
    rotate_x += 15.0 * (sin(p) - sin(q));

    true_redisplay = 1;
  }
  
  if(true_redisplay) {
    glLoadIdentity();
 
   /* Define viewing transformation */
    gluLookAt((GLdouble)eye[0],(GLdouble)eye[1],(GLdouble)eye[2],
	      (GLdouble)center[0],(GLdouble)center[1],(GLdouble)center[2],
	      (GLdouble)up[0],(GLdouble)up[1],(GLdouble)up[2]);

    glTranslatef (center[0], center[1], center[2] + translate_z);
    glRotatef(rotate_x, 1.0, 0.0, 0.0);
    glRotatef(rotate_y, 0.0, 1.0, 0.0);
    glTranslatef(-center[0], -center[1], -center[2]);

    true_redisplay = 0;
  }

  render_sand_surface();

  for(unsigned i = 0; i < 2 ; i++) {
    glBegin (GL_LINE_LOOP);
    glColor4f (1.0f, 1.0f, 1.0f, 0.9f);
    glVertex3f (min_ext[0], min_ext[1], i ? max_ext[2] : min_ext[2]);
    glVertex3f (max_ext[0], min_ext[1], i ? max_ext[2] : min_ext[2]);
    glVertex3f (max_ext[0], max_ext[1], i ? max_ext[2] : min_ext[2]);
    glVertex3f (min_ext[0], max_ext[1], i ? max_ext[2] : min_ext[2]);
    glEnd ();
  }

#ifndef PIX_SURFACE
  if(texture) {
    glEnable(GL_TEXTURE_2D);     // Enable 2D texturing

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, SAND_DIM, SAND_DIM, 
		 0, GL_RGB , GL_FLOAT, texture);

    glBegin(GL_QUADS);
    glTexCoord2i(0, 0);
    glVertex3f(min_ext[0], max_ext[1], min_ext[2]);glNormal3f(0.0, 1.0, 0.0);
    glTexCoord2i(0, 1);
    glVertex3f(min_ext[0], max_ext[1], max_ext[2]);glNormal3f(0.0, 1.0, 0.0);
    glTexCoord2i(1, 1);
    glVertex3f(max_ext[0], max_ext[1], max_ext[2]);glNormal3f(0.0, 1.0, 0.0);
    glTexCoord2i(1, 0);
    glVertex3f(max_ext[0], max_ext[1], min_ext[2]);glNormal3f(0.0, 1.0, 0.0);
    glEnd();
    glDisable(GL_TEXTURE_2D);     // Enable 2D texturing
  }
#endif
}

void do_display ()
{
  glClear (GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  display_draw_scene ();
  glutSwapBuffers ();
}

extern void sable_animate (unsigned iterations);

volatile int keepCool = 0;


#define TIME_DIFF(t1, t2) \
    ((t2.tv_sec - t1.tv_sec) * 1000000 + (t2.tv_usec - t1.tv_usec))



volatile static int nbFrames = 0;
volatile int computeTime=0;

int periods[]={2000000,1000000,500000,250000,100000,50000,10000};
#define MaxDisplayPeriod (sizeof(periods)/sizeof(periods[0]))
volatile int displayPeriod = MaxDisplayPeriod ; // 10 fps

volatile int nbIterations = 0;
int iterations[]={1,2,5,10,25,50,100,250,500,1000,10000};
#define MaxNbIterations (sizeof(iterations)/sizeof(iterations[0]))


void printFPS()
{
  static struct timeval lastTime;
  static int jamais_init = 1;
  struct timeval now;

  if (jamais_init)
    {
      jamais_init = 0;
      gettimeofday(&lastTime,NULL);
    }

  gettimeofday(&now,NULL);
  
  if ( nbFrames > 0 && TIME_DIFF(lastTime,now) >= 1000000 ){ // If last prinf() was more than 1 sec ago
    // printf and reset timer
    printf("\r %.1f ms/frame Compute Time %.3f ms/frame ",
	   (TIME_DIFF(lastTime,now))/(double)(nbFrames*1000),computeTime/(double)(nbFrames*1000) );
    printf("%d iterations à la fois", iterations[nbIterations]);
    fflush(stdout);
    nbFrames = 0;
    computeTime = 0;
    lastTime = now;
  }
}

void updateDisplay(int i)
{
  //  printf("ouf!\n");
  keepCool = 0;
  gettimeofday (&lastDisplayTimeval,NULL);
  glutPostRedisplay();
}

void idle(void)
{
  float *colors;
  struct timeval now;
  static int firstCall = 1;
  
  printFPS();

  if(!keepCool)
    {
      struct timeval t1,t2;
      gettimeofday (&t1,NULL);
      colors = compute_func_ptr (iterations[nbIterations]);
      gettimeofday (&t2,NULL);
      computeTime += TIME_DIFF(t1,t2);
      nbFrames++;

#ifdef PIX_SURFACE
      if (colors && colors != DYNAMIC_COLORING)
	texture = colors;

      sand_surface_refresh (colors != STATIC_COLORING);

      if (colors) {
	// Refresh colors
	glBindBuffer(GL_ARRAY_BUFFER, s_vbocid);
	glBufferSubData(GL_ARRAY_BUFFER, 0, s_vertices*3*sizeof(float),
			s_vbo_color);
      }
#else
      sand_surface_refresh (colors == DYNAMIC_COLORING);
      
      if (colors) {
	// Refresh colors
	glBindBuffer(GL_ARRAY_BUFFER, s_vbocid);
	glBufferSubData(GL_ARRAY_BUFFER, 0, s_vertices*3*sizeof(float),
			(colors == DYNAMIC_COLORING) ? s_vbo_color : colors);
	
	texture = (colors == DYNAMIC_COLORING) ? s_vbo_color : colors;
      }
#endif
      
      // Refresh vertices
      glBindBuffer(GL_ARRAY_BUFFER, s_vbovid);
      glBufferSubData(GL_ARRAY_BUFFER, 0, s_vertices*3*sizeof(float), s_vbo_vertex);
      
      
      gettimeofday (&now,NULL);
      int timediff = TIME_DIFF(lastDisplayTimeval,now);
      if (timediff>=(periods[displayPeriod]) || firstCall )
	{
	  firstCall = 0;
	  lastDisplayTimeval = now; 
	  glutPostRedisplay ();
	}
      else
	{
	  keepCool = 1;
	  //	  printf("displayPeriod-timediff  %d \n",(periods[displayPeriod]-timediff)/1000); 
	  glutTimerFunc((1000+periods[displayPeriod]-timediff)/1000, updateDisplay,0);
	}
    }
}
  void initView (float *min_ext, float *max_ext)
  {
    GLfloat light_diffuse[]   = {1.0, 1.0, 1.0, 1.0};
  GLfloat light_position[] = {0.5, 0.5, 1.0, 0.0};
  float dif_ext[3];
  int i;

  /* Define normal light */
  glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
  glLightfv(GL_LIGHT0, GL_POSITION, light_position);

  /* Enable a single OpenGL light */
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);

  /* Use depth buffering for hidden surface elimination */
  glEnable(GL_DEPTH_TEST);

  glEnable(GL_POINT_SPRITE);
  glEnable(GL_POINT_SMOOTH);
  glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);

  /* get diagonal and average distance of extent */
  for (i=0; i<3; i++)
    dif_ext[i] = max_ext[i] - min_ext[i];
  
  dis = 0.0;
  for (i = 0; i < 3; i++)
    dis = MAX(dis, dif_ext[i]);
  dis *= 1.2;

  /* set center in world space */
  for (i=0; i<3; i++)
    center[i] = min_ext[i] + dif_ext[i] * 0.5;

  /* set initial eye & look at location in world space */
  eye[0] = center[0];
  eye[1] = center[1];
  eye[2] = center[2] + dis;
  up[0] = 0.0;
  up[1] = 1.0;
  up[2] = 0.0;

  // set parameters for gluPerspective()
  near = dis-0.5*dif_ext[2];
  near_clip = (GLdouble) (0.5 * near);
  far_clip = (GLdouble) (2.0 * (dis + 0.5 * dif_ext[2]));

  // Field of view 
  fovy = (GLdouble) (0.5 * dif_ext[1] / near);
  fovy = (GLdouble) (2 * atan ((double) fovy) / M_PI * 180.0);

  // Enable the color material mode
  glEnable(GL_COLOR_MATERIAL);
}

void appDestroy()
{
  //atomFinalize();
  vboFinalize();

  if(glutWindowHandle)
    glutDestroyWindow(glutWindowHandle);

  exit(0);
}


void appKeyboard(unsigned char key, int x, int y)
{
    //this way we can exit the program cleanly
    switch(key)
    {
    case '<' : translate_z -= 0.1; true_redisplay = 1; glutPostRedisplay(); break;
    case '>' : translate_z += 0.1; true_redisplay = 1; glutPostRedisplay(); break;
    case 'r':
    case 'R': camera_rotate = 1 - camera_rotate; break;
    case 's':
    case 'S': SAND_MAX_HEIGHT = max_value; printf ("Setting MAX_HEIGHT to %d\n", max_value); break;
#ifdef PIX_SURFACE
    case 'f':
    case 'F': flatten_surface = 1 - flatten_surface; break;
#endif
    case '\033': // escape quits
    case '\015': // Enter quits    
    case 'Q':    // Q quits
    case 'q':    // q (or escape) quits
      // Cleanup up and quit
      appDestroy();     break;
    }
}


void appSpecialKeyboard(int key, int x, int y)
{
  //this way we can exit the program cleanly
  switch(key)
    {
    case GLUT_KEY_UP:
      if (displayPeriod < MaxDisplayPeriod)
	displayPeriod++;
	  printf ("Setting display period to %d \n", displayPeriod);
	  break;
    case GLUT_KEY_DOWN:
      if (displayPeriod > 0)
	displayPeriod--;
      printf ("Setting display period to %d \n", displayPeriod);
      break;
    case GLUT_KEY_RIGHT:
      if (nbIterations < MaxNbIterations)
	nbIterations++;
      printf ("Setting iteration to %d \n", nbIterations);
      break;
    case GLUT_KEY_LEFT:
      if (nbIterations > 0)
	nbIterations--;
      printf ("Setting iteration to %d \n", nbIterations);
      break;
    }
}



void appMouse(int button, int state, int x, int y)
{
    //handle mouse interaction for rotating/zooming the view
    if (state == GLUT_DOWN) {
        mouse_buttons |= 1<<button;
    } else if (state == GLUT_UP) {
        mouse_buttons = 0;
    }

    mouse_old_x = x;
    mouse_old_y = y;
}

void appMotion(int x, int y)
{
  float dx, dy;

  dx = x - mouse_old_x;
  dy = y - mouse_old_y;

  if (mouse_buttons & 1) {
    rotate_x += dy * 0.2;
    rotate_y += dx * 0.2;
  } else if (mouse_buttons & 4) {
    translate_z += dy * 0.1;
  }

  mouse_old_x = x;
  mouse_old_y = y;

  true_redisplay = 1;
  glutPostRedisplay();
}

void display_init (int argc, char **argv, unsigned dim, unsigned max_height,
		   get_func_t get_func,
		   compute_func_t compute_func)
{
  SAND_DIM = dim;
  SAND_MAX_HEIGHT = max_height;
  get_func_ptr = get_func;
  compute_func_ptr = compute_func;
  
  glutInit(&argc, argv);

  /* Initialize arrays of vertices and triangles */
  initializeVBO();

  /* Set up an window */
  /* Initialize display mode */
  glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH);
  /* Specify window size */
  glutInitWindowSize(winx, winy);
  /* Open window */
  glutWindowHandle = glutCreateWindow("Tas de sable");

  /* Initialize view */
  initView(min_ext, max_ext);

  /* Set a glut callback functions */
  glutDisplayFunc(do_display);
  glutReshapeFunc(reshape);
  glutIdleFunc(idle);

  glutKeyboardFunc(appKeyboard);
  glutSpecialFunc(appSpecialKeyboard);
 glutMouseFunc(appMouse);
  glutMotionFunc(appMotion);

  buildVBO();

  /* Start main display loop */
  glutMainLoop();
}
