#include <iostream>

#include "game.hpp"
#include <GL/glut.h>
#include <thread>
#include <chrono>

    using namespace std::chrono;

    milliseconds last_frame_time;

game G;

  int wd;                   /* GLUT window handle */
  GLdouble width, height;   /* window width and height */

  bool flag_display;

using std::cout;
using std::endl;

int frame_number;

void init()
{
  width  = 500.0;                 /* initial window width and height, */
  height = 500.0;                  /* within which we draw. */
}

/* Callback functions for GLUT */

/* Draw the window - this is where all the GL actions are */
void display(void)
{
  //if(flag_display) return;
  //flag_display=true;

  const double fps=25.0;
  milliseconds current_time = duration_cast<milliseconds>(high_resolution_clock::now().time_since_epoch());
  milliseconds time_past_from_last_frame = current_time-last_frame_time;

  //if(time_past_from_last_frame<milliseconds(100)) return;

  milliseconds frame_duration(int(1000.0/fps));
  milliseconds time_residual=frame_duration-time_past_from_last_frame;

  //cout<<"sleep for "<<time_residual.count()<<endl;
  std::this_thread::sleep_for(time_residual);
  //std::this_thread::sleep_for(frame_duration);


  last_frame_time = duration_cast<milliseconds>(high_resolution_clock::now().time_since_epoch());
  //cout<<"Frame "<<frame_number++<<endl;

  /*for(int k=0; k<80; k++)
  {
    G.PS->step();
  }

  G.PS->INT->status(cout);*/

  /* clear the screen to white */
  glClear(GL_COLOR_BUFFER_BIT);

  G.R->draw_particles();
  G.R->draw_frame();

  glFlush();

  glutSwapBuffers();


  //flag_display=false;
}


/* Called when window is resized,
   also when window is first created,
   before the first call to display(). */
void
reshape(int w, int h)
{
  /* save new screen dimensions */
  width = (GLdouble) w;
  height = (GLdouble) h;

  /* tell OpenGL to use the whole window for drawing */
  glViewport(0, 0, (GLsizei) width, (GLsizei) height);

  /* do an orthographic parallel projection with the coordinate
     system set to first quadrant, limited by screen/window size */
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0.0, width, 0.0, height, -1.f, 1.f);

  return;
}

void
kbd(unsigned char key, int x, int y)
{
  switch((char)key) {
  case 'q':
  case 27:    /* ESC */
    glutDestroyWindow(wd);
    exit(0);
  case 'u':
    display();
  default:
    break;
  }

  return;
}

void mouse_move(int x, int y) {
  vect c;
  c.x = -1. + 2. * (x / width);
  c.y = 1. - 2. * (y / height);

  G.PS->SetForce(c);
}

void mouse(int button, int state, int x, int y)
{
  cout<<button<<" "<<state<<" "<<x<<" "<<y<<endl;
  G.PS->SetForce(state == 0);
}


void cycle()
{
  cout<<"Computation thread started"<<endl;
  while(true)
  {
    G.PS->step();
  }
}



int main(int argc, char *argv[])
{
    frame_number=0;
    /* perform initialization NOT OpenGL/GLUT dependent,
       as we haven't created a GLUT window yet */
    init();

    /* initialize GLUT, let it extract command-line
       GLUT options that you may provide
       - NOTE THE '&' BEFORE argc */
    glutInit(&argc, argv);

    /* specify the display to be single
       buffered and color as RGBA values */
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA | GLUT_DOUBLE );

    /* set the initial window size */
    glutInitWindowSize((int) width, (int) height);

    /* create the window and store the handle to it */
    wd = glutCreateWindow("part-toy" /* title */ );

    /* --- register callbacks with GLUT --- */

    /* register function to handle window resizes */
    glutReshapeFunc(reshape);

    /* register keyboard event processing function */
    glutKeyboardFunc(kbd);

    glutMouseFunc(mouse);
    glutMotionFunc(mouse_move);

    /* register function that draws in the window */
    glutDisplayFunc(display);
    glutIdleFunc(glutPostRedisplay);

    /* init GL */
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glColor3f(0.0, 0.0, 0.0);
    glLineWidth(3.0);

    std::thread computation_thread(cycle);

    flag_display=false;

    /* start the GLUT main loop */
   glutMainLoop();

}
