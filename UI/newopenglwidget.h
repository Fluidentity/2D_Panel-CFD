#ifndef NEWOPENGLWIDGET_H
#define NEWOPENGLWIDGET_H

#include <QOpenGLWidget>
#include <QTimer>


class NewOpenGLWidget : public QOpenGLWidget
{
public:
  explicit NewOpenGLWidget(QWidget *parent = 0);

  void initializeGL();
  void paintGL();
  void resizeGL(int w, int h);

  float D;
  float x=0, y=0, z=0;
  bool wired;
  float ang;
private:
  QTimer timer;
};

#endif // NEWOPENGLWIDGET_H
