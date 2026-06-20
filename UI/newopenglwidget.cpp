#include "newopenglwidget.h"
#include <GL/glut.h>

NewOpenGLWidget::NewOpenGLWidget(QWidget *parent)
{
  this->D = 1;
  this->wired = true;
  this->ang = 0.5;

  connect(&timer, SIGNAL(timeout()), this, SLOT(update()));
  timer.start(16);

}

void NewOpenGLWidget::initializeGL()
{
  glClearColor(0.3, 0.3, 0.3, 1);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHTING);
  glEnable(GL_COLOR_MATERIAL);

}

void NewOpenGLWidget::paintGL()
{
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt( 0,0,5, 0,0,0, 0,1,0 );
  ang +=0.5;

  glRotatef(ang, 1, 1, 1);
  glTranslatef(x, y, z);

  glColor3f(1, 0, 0);
  if(wired)
    glutWireSphere(D,20,20);
  else
    glutSolidSphere(D,20,20);


}

void NewOpenGLWidget::resizeGL(int w, int h)
{
  glViewport(0, 0, w, h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(45.0, (float)w/h, 0.01, 100.0);
  update();


}
