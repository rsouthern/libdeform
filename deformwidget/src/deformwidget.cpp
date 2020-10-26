#include "deformwidget.h"

// Include for deform package
#include <laplacian_COT.h>
#include <laplacian_DUAL.h>
#include <laplacian_DUAL2.h>
#include <cholmod_shared.h>


DeformWidget::DeformWidget(QWidget* parent) : DeformWidgetParent(parent) {
    installEventFilter(this);
    currentDeformer = 0;

    unique_vertex_colors = NULL;

    add_draw_mode("Vertex Selection");
    add_draw_mode("Front Face Vertex Selection");

    selectionFrontFace = true;
    fboInit = false;
    state = STATE_IDLE;
    op = NULL;
}


DeformWidget::~DeformWidget() {
    if (unique_vertex_colors != NULL) free(unique_vertex_colors);
}

void DeformWidget::open_mesh_gui(QString fname) {
    OpenMesh::Utils::Timer t;
    t.start();
    if ( fname.isEmpty() || !open_mesh(fname.toLocal8Bit(), _options) )
    {
        QString msg = "Cannot read mesh from file:\n '";
        msg += fname;
        msg += "'";
        QMessageBox::critical( NULL, windowTitle(), msg);
    }
    t.stop();
    std::cout << "Loaded mesh in ~" << t.as_string() << std::endl;
}

void DeformWidget::open_texture_gui(QString fname) {
    if ( fname.isEmpty() || !open_texture( fname.toLocal8Bit() ) )
    {
        QString msg = "Cannot load texture image from file:\n '";
        msg += fname;
        msg += "'\n\nPossible reasons:\n";
        msg += "- Mesh file didn't provide texture coordinates\n";
        msg += "- Texture file does not exist\n";
        msg += "- Texture file is not accessible.\n";
        QMessageBox::warning( NULL, windowTitle(), msg );
    }
}

void DeformWidget::query_open_mesh_file() {
    QString fileName = QFileDialog::getOpenFileName(this,
                                                    tr("Open mesh file"),
                                                    tr(""),
                                                    tr("OBJ Files (*.obj);;"
                                                       "OFF Files (*.off);;"
                                                       "STL Files (*.stl);;"
                                                       "All Files (*)"));
    if (!fileName.isEmpty())
        open_mesh_gui(fileName);
}

void DeformWidget::query_open_texture_file() {
    QString fileName = QFileDialog::getOpenFileName(this,
                                                    tr("Open texture file"),
                                                    tr(""),
                                                    tr("PNG Files (*.png);;"
                                                       "BMP Files (*.bmp);;"
                                                       "GIF Files (*.gif);;"
                                                       "JPEG Files (*.jpg);;"
                                                       "TIFF Files (*.tif);;"
                                                       "All Files (*)"));
    if (!fileName.isEmpty())
        open_texture_gui(fileName);
}

void DeformWidget::set_front_face(bool isOn) {
    selectionFrontFace = isOn;
}

void DeformWidget::set_deformer() {
    QDialog *dDialog = new QDialog(this);
    dDialog->setWindowTitle("Set Deformer Type");
    QButtonGroup *dButtons = new QButtonGroup(dDialog);
    QRadioButton *l = new QRadioButton(dDialog->tr("&Equal Weight Laplacian"), dDialog);
    QRadioButton *l_cot = new QRadioButton(dDialog->tr("&Cotangent Weight Laplacian"), dDialog);
    QRadioButton *l_dual = new QRadioButton(dDialog->tr("&Dual Laplacian"), dDialog);
    QRadioButton *l_dual2 = new QRadioButton(dDialog->tr("&Second Order Dual Laplacian"), dDialog);

    QVBoxLayout *layout = new QVBoxLayout;
    dButtons->addButton(l,DEF_LAPLACIAN);
    dButtons->addButton(l_cot,DEF_LAPLACIAN_COT);
    dButtons->addButton(l_dual,DEF_LAPLACIAN_DUAL);
    dButtons->addButton(l_dual2,DEF_LAPLACIAN_DUAL2);

    layout->addWidget(l);
    layout->addWidget(l_cot);
    layout->addWidget(l_dual);
    layout->addWidget(l_dual2);

    switch(currentDeformer) {
    case DEF_LAPLACIAN_COT:
        l_cot->setChecked(true);
        l_cot->setFocus();
        break;
    case DEF_LAPLACIAN_DUAL:
        l_dual->setChecked(true);
        l_dual->setFocus();
        break;
    case DEF_LAPLACIAN_DUAL2:
        l_dual2->setChecked(true);
        l_dual2->setFocus();
        break;
    default: // DEF_LAPLACIAN
        l->setChecked(true);
        l->setFocus();
        break;
    }

    // Add a button at the bottom to accept
    QPushButton *dAccept = new QPushButton("Ok", dDialog);
    QObject::connect(dAccept,SIGNAL(clicked()),dDialog,SLOT(accept()));
    layout->addWidget(dAccept);
    dDialog->setLayout(layout);


    if (dDialog->exec() == QDialog::Accepted) {
        fprintf(stderr, "\nCurrently selected method = %d", dButtons->checkedId());
        currentDeformer = dButtons->checkedId();
    }
}

/**
  * Clear away the list of constraints and anchors
  */
void DeformWidget::clear_deformer() {
    anchors.clear();
    constraints.clear();
}


/**
  * Our own open_mesh() function. This updates the CHOLMOD data stored in the deformer widget.
  */
bool DeformWidget::open_mesh(const char* _filename, IO::Options _opt) {
    bool ret_val = MeshViewerWidgetT<Mesh>::open_mesh(_filename,_opt);

    if (ret_val) {
        // Added this here to update the Cholmod data
        mesh_.updateX();
        mesh_.updateTri();
        anchors.resize(mesh_.n_vertices(),false);
        constraints.resize(mesh_.n_vertices(),false);
        generate_unique_vertex_colors();
    }
    return ret_val;
}

QPointF DeformWidget::pixelPosToViewPos(const QPointF& p) {
    return QPointF(2.0 * float(p.x()) / width() - 1.0,
                   1.0 - 2.0 * float(p.y()) / height());
}



bool DeformWidget::eventFilter(QObject *o, QEvent *e) {
    QKeyEvent *ke;
    QMouseEvent *me;
    QWheelEvent *we;
    int i,j;

    GLfloat minCorner[3] = {0.0f,0.0f,0.0f};
    GLfloat maxCorner[3] = {0.0f,0.0f,0.0f};
    GLfloat midPt[3] = {0.0f,0.0f,0.0f};

    switch(e->type()) {
    case (QEvent::FocusIn):
        // When we focus into the object we need to retrieve the key state
        // to ensure that we are not deforming or specifying anchors
        break;
    case (QEvent::FocusOut):
        // When we leave focus, we need to cancel any existing deforming /
        // selection states.
        state = STATE_IDLE;
        sbox.type = STYPE_NOTHING;
        break;
    case (QEvent::MouseButtonPress):
        me = (QMouseEvent*) e;
        switch(me->button()) {
        case (Qt::LeftButton):
            deformTrackBall.MouseButton(0,1,me->x(),me->y());
            break;
        case (Qt::MidButton):
            deformTrackBall.MouseButton(1,1,me->x(),me->y());
            break;
        case (Qt::RightButton):
            deformTrackBall.MouseButton(2,1,me->x(),me->y());
            break;
        default:
            break;
        }
        break;
    case (QEvent::MouseButtonRelease):
        me = (QMouseEvent*) e;
        switch(me->button()) {
        case (Qt::LeftButton):
            deformTrackBall.MouseButton(0,0,me->x(),me->y());
            break;
        case (Qt::MidButton):
            deformTrackBall.MouseButton(1,0,me->x(),me->y());
            break;
        case (Qt::RightButton):
            deformTrackBall.MouseButton(2,0,me->x(),me->y());
            break;
        default:
            break;
        }
        break;
    case (QEvent::Wheel):
        we = (QWheelEvent*) e;
        if (state == STATE_SELECT) {
            sbox.radius = max(1, sbox.radius + (we->delta()/60));
            updateGL();
            return true;
        }
        break;
    case (QEvent::MouseMove):
        // This does nothing unless we are deforming or selecting anchors
        me = (QMouseEvent*) e;
        switch(state) {
        case(STATE_SELECT):
            // If the left button is down we are selecting anchors
            if (me->buttons() & Qt::LeftButton) {
                fprintf(stderr,"\nSTATE_SELECT Anchors");                
                sbox.x = me->x(); sbox.y = height() - me->y();
                sbox.type = STYPE_ANCHORS;
                find_selected_vertices(me->x(), height() - me->y(), sbox.radius, &anchors);
            } else if (me->buttons() & Qt::MidButton) {
                fprintf(stderr,"\nSTATE_SELECT Constraints");
                sbox.x = me->x(); sbox.y = height() - me->y();
                sbox.type = STYPE_CONSTRAINTS;
                find_selected_vertices(me->x(), height() - me->y(), sbox.radius, &constraints);
            } else {
                sbox.type = STYPE_NOTHING;
            }
            return true;
            break;
        case(STATE_DEFORM):
            fprintf(stderr,"\nSTATE_DEFORM");

            if (op != NULL) {
                // Create a transformation matrix for deformation
                deformTrackBall.MouseMove(me->x(), me->y());
                op->transformHandles(deformTrackBall.transformMat());
                op->update_cuda();

                // Update the data in a mesh file for the purposes of saving
                mesh_.setX(op->getX());

                deformTrackBall.reset();

                // Redraw the screen
                updateGL();
                return true;
            }
            break;
        default: // Idle
            break;
        }
        break;
    case (QEvent::KeyPress):
        // Determine if we are changing our deformation status
        ke = (QKeyEvent *) e;
        switch(ke->key()) {
        case KEY_SELECT:
            if ((state == STATE_IDLE) && (!ke->isAutoRepeat())) state = STATE_SELECT;
            break;
        case KEY_DEFORM:
            if ((state == STATE_IDLE) && (!ke->isAutoRepeat())) {
                printf("\nEntering Deformation Mode");
                // Initialise a deformer with our currenty selected constraints
                if (op != NULL) delete op;
                switch(currentDeformer) {
                case (DEF_LAPLACIAN_COT):
                    op = new Laplacian_COT();
                    break;
                case (DEF_LAPLACIAN_DUAL):
                    op = new Laplacian_DUAL();
                    break;
                case (DEF_LAPLACIAN_DUAL2):
                    op = new Laplacian_DUAL2();
                    break;
                default:
                    op = new Laplacian();
                    break;
                }
                if (op != NULL) {
                    for (i=0; i<mesh_.n_vertices(); ++i) {
                        // append anchors and constraints
                        if (anchors[i])
                            op->addAnchor(i);
                        if (constraints[i])
                            op->addHandle(i);

                        // Determine the midpoint of the handles
                        if (constraints[i]) {
                            for (j=0; j<3; ++j) {
                                if (minCorner[j] > mesh_.points()[i][j]) minCorner[j] = mesh_.points()[i][j];
                                if (   maxCorner[j] < mesh_.points()[i][j]) maxCorner[j] = mesh_.points()[i][j];
                            }
                        }
                    }

                    FILE *fid = fopen("constraints.dat", "w");
                    if (fid) {
                        op->saveConstraints(fid);
                        fclose(fid);
                    }
                    op->initialise(&mesh_);

                    // Set the origin of our deformation trackball
                    midPt[0] = 0.5f*(maxCorner[0]-minCorner[0]);
                    midPt[1] = 0.5f*(maxCorner[1]-minCorner[1]);
                    midPt[2] = 0.5f*(maxCorner[2]-minCorner[2]);
                    deformTrackBall.setModelOrigin(midPt[0], midPt[1], midPt[2]);
                    deformTrackBall.setModelDistance(maxCorner[0] - midPt[0]); // Just use x manhattan distance
                    deformTrackBall.setWindowDimensions(width(), height());
                    state = STATE_DEFORM;
                }
            }

            break;
        default:
            break;
        }
        break;
    case (QEvent::KeyRelease):
        // Determine if we are changing our deformation status
        ke = (QKeyEvent *) e;
        switch(ke->key()) {
        case KEY_SELECT:
            if ((state == STATE_SELECT) && !ke->isAutoRepeat()) state = STATE_IDLE;
            break;
        case KEY_DEFORM:
            if ((state == STATE_DEFORM) && !ke->isAutoRepeat()) {
                state = STATE_IDLE;
                printf("\nLeaving Deformation Mode");
            }
            break;
        default:
            break;
        }
        break;
    default:
        break;
    }
    return DeformWidgetParent::eventFilter(o,e);
}

/**
  * If we were organised this would be a trait of the OpenMesh structure.
  */
void DeformWidget::generate_unique_vertex_colors() {
    Mesh::VertexIter v_it;

    if (unique_vertex_colors != NULL) free(unique_vertex_colors);
    unique_vertex_colors = (unsigned char*) malloc(sizeof(unsigned char)*3*mesh_.n_vertices());

    GLuint i, idx;

    for (v_it = mesh_.vertices_begin(),i=0; v_it != mesh_.vertices_end(); ++v_it,++i) {
        // Map the index to a 3f color value
        cmap.idxToColor(i, &(unique_vertex_colors[3*i+0]),
                            &(unique_vertex_colors[3*i+1]),
                            &(unique_vertex_colors[3*i+2]));

        // Sanity check
        cmap.colorToIdx(unique_vertex_colors[3*i+0], unique_vertex_colors[3*i+1], unique_vertex_colors[3*i+2], &idx);

        if (idx != i) {
            fprintf(stderr, "\nDeformWidget::generate_unique_vertex_colors() - cmap[%u,%u,%u] maps to %u not %u!",
                    unique_vertex_colors[3*i+0], unique_vertex_colors[3*i+1], unique_vertex_colors[3*i+2], idx, i);
            exit(0);
        }
    }
}

/**
  * Need to override the draw_scene function so we can catch the vertex selection action
  */
void DeformWidget::draw_scene(const std::string& _draw_mode) {
    // Make sure there is something there to select
    if (mesh_.n_vertices() == 0) return;

    // Initialise an FBO here (we can't override the initializeGL without changing the parent class!)
    if (!fbo.isInit()) fbo.init(width(), height());

    // Render the selection into the frame buffer
    fbo.bind();
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glDisable(GL_LIGHTING);
    glDisable (GL_BLEND);
    glDisable (GL_DITHER);
    glDisable (GL_FOG);
    glDisable (GL_TEXTURE_1D);
    glDisable (GL_TEXTURE_2D);
    glDisable (GL_TEXTURE_3D);
    glShadeModel (GL_FLAT);
    if (selectionFrontFace) {
        // For front face selection draw black polygons with an offset and then
        // the colored vertices
        glEnable(GL_POLYGON_OFFSET_FILL);
        glPolygonOffset(2.0f,1.0f);
        draw_openmesh("Filled Black Polygons");
        draw_openmesh("Unique Vertex Colors");
        glDisable(GL_POLYGON_OFFSET_FILL);
    } else {
        // Just draw the colored vertices without anything else
        draw_openmesh("Unique Vertex Colors");
    }
    // Unbind the frame buffer
    fbo.unbind();

    glClearColor(1.0, 1.0, 1.0, 1.0);

    if (_draw_mode == "Vertex Selection") {
        glDisable(GL_LIGHTING);
        draw_openmesh("Unique Vertex Colors");
    } else if (_draw_mode == "Front Face Vertex Selection") {
        glDisable(GL_LIGHTING);
        glEnable(GL_POLYGON_OFFSET_FILL);
        glPolygonOffset(2.0f,1.0f);
        draw_openmesh("Filled Black Polygons");
        draw_openmesh("Unique Vertex Colors");
        glDisable(GL_POLYGON_OFFSET_FILL);
    }

    // Draw the trackball
    if (state == STATE_DEFORM)
        deformTrackBall.draw();


    // We will now always display our constraints and anchors
    glEnable(GL_LIGHTING);
    glEnable(GL_BLEND);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(2.0f,1.0f);
    DeformWidgetParent::draw_scene(_draw_mode);

    glDisable(GL_LIGHTING);
    draw_openmesh("Anchor Vertices");
    draw_openmesh("Constraint Vertices");
    draw_openmesh("Selection Box");
    glEnable(GL_LIGHTING);
    glDisable(GL_POLYGON_OFFSET_FILL);
}

/**
  * Render an image to determine mouse selection.
  */
void DeformWidget::draw_openmesh(const std::string& _draw_mode) {
    //fprintf(stderr,"\nDeformWidget::draw_openmesh(%s)", _draw_mode.c_str());
    if (_draw_mode == "Filled Black Polygons") {
        glColor3f(0.0f,0.0f,0.0f);
        glEnableClientState(GL_VERTEX_ARRAY);
        glVertexPointer(3, GL_FLOAT, 0, mesh_.points());

        glEnableClientState(GL_NORMAL_ARRAY);
        glNormalPointer(GL_FLOAT, 0, mesh_.vertex_normals());

        MyStripifier::StripsIterator strip_it = strips_.begin();
        MyStripifier::StripsIterator strip_last = strips_.end();

        // Draw all strips
        for (; strip_it!=strip_last; ++strip_it) {
            glDrawElements(GL_TRIANGLE_STRIP,
                           strip_it->size(), GL_UNSIGNED_INT, &(*strip_it)[0] );
        }
        glDisableClientState(GL_VERTEX_ARRAY);
        glDisableClientState(GL_NORMAL_ARRAY);
        glDisableClientState(GL_TEXTURE_COORD_ARRAY);
    } else if (_draw_mode == "Unique Vertex Colors") {
        // Might need these to ensure the things are visible
        glPointSize(1.0f);

        // Render the vertices into the current buffer
        glEnableClientState(GL_VERTEX_ARRAY);
        glVertexPointer(3, GL_FLOAT, 0, mesh_.points());

        // Use our unique vertex colors
        glEnableClientState(GL_COLOR_ARRAY);
        glColorPointer(3, GL_UNSIGNED_BYTE, 0, unique_vertex_colors);

        // Draw all of it
        glDrawArrays( GL_POINTS, 0, mesh_.n_vertices() );

        // Disable the array stuff
        glDisableClientState(GL_VERTEX_ARRAY);
        glDisableClientState(GL_COLOR_ARRAY);
        glPointSize(1.0f);
    } else if (_draw_mode == "Anchor Vertices") {
        glColor3f(0.0f,0.0f,1.0f);
        glPointSize(4.0f);
        glBegin(GL_POINTS);
        uint i;
        for (i=0; i<anchors.size(); ++i) {
            if (anchors[i]) {
                glVertex3f(mesh_.points()[i][0], mesh_.points()[i][1], mesh_.points()[i][2]);
            }
        }
        glEnd(); //GL_POINTS
        glPointSize(1.0f);

    } else if (_draw_mode == "Constraint Vertices") {
        glColor3f(1.0f,0.0f,0.0f);
        glPointSize(4.0f);
        list<uint>::iterator it;
        glBegin(GL_POINTS);
        uint i;
        for (i=0; i<constraints.size(); ++i) {
            if (constraints[i]) {
                glVertex3f(mesh_.points()[i][0], mesh_.points()[i][1], mesh_.points()[i][2]);
            }
        }
        glEnd(); //GL_POINTS
        glPointSize(1.0f);
    } else if (_draw_mode == "Selection Box") {
        // Make a 2D projection
        glMatrixMode(GL_PROJECTION);
        glPushMatrix();
        glLoadIdentity();
        glOrtho( -GLdouble(width())*0.5 , GLdouble(width())*0.5,
                 -GLdouble(height())*0.5, GLdouble(height())*0.5,
                 -10.0f, 10.0f );

        // Put the transform here to position the box in the centre
        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        glLoadIdentity();
        glTranslatef(-GLfloat(width())*0.5, -GLfloat(height())*0.5, 0.0f);

        // Draw something
        sbox.draw();

        glPopMatrix();

        // Revert to original projection
        glMatrixMode(GL_PROJECTION);
        glPopMatrix();
        glMatrixMode(GL_MODELVIEW);
    } else {
        return DeformWidgetParent::draw_openmesh(_draw_mode);
    }
}

void SelectionBox::draw() {
    switch(type) {
    case (STYPE_ANCHORS):
        glColor4f(0.0f, 1.0f, 0.0f, 0.5f);
        break;
    case (STYPE_CONSTRAINTS):
        glColor4f(0.0f, 0.0f, 1.0f, 0.5f);
        break;
    default:
        return;
    }
    GLfloat minx = GLfloat(x) - GLfloat(radius);
    GLfloat maxx = GLfloat(x) + GLfloat(radius);
    GLfloat miny = GLfloat(y) - GLfloat(radius);
    GLfloat maxy = GLfloat(y) + GLfloat(radius);

    glBegin(GL_LINE_LOOP); {
        glVertex3f(minx, miny, 0.0f);
        glVertex3f(minx, maxy, 0.0f);
        glVertex3f(maxx, maxy, 0.0f);
        glVertex3f(maxx, miny, 0.0f);
    } glEnd();
}

/**
  * Retrieve any vertices under the selection region and set it in the passed list. If
  * there are none, reset all vertices in the list. This will use glReadPixels from the
  * existing selection buffer, scan the pixels in the region and use the color map to
  * retrieve the selected vertices.
  */
void DeformWidget::find_selected_vertices(uint x, uint y, uint radius, vector<bool> *vvec) {
    // Determine the extents of the selection region based on the current resolution
    uint minx, miny;
    uint maxx = x+radius;
    uint maxy = y+radius;
    if (x < radius) {
        minx = 0;
    } else {
        minx = x-radius;
    }
    if (y < radius) {
        miny = 0;
    } else {
        miny = y-radius;
    }
    if (maxx > fbo.width()) maxx = fbo.width()-1;
    if (maxy > fbo.height()) maxy = fbo.height()-1;

    // Size of image data = 2r * 2r * 3 (color channels)

    GLuint w = maxx - minx;
    GLuint h = maxy - miny;

    // A zero selection size isn't particularly interesting
    if (w == 0 || h == 0) return;

    unsigned char *imageData = new unsigned char[3*w*h];

    // Get the chunk of memory from the graphics card    
    makeCurrent();

    if (fbo.bind()) {
        glReadPixels(minx, miny, w, h, GL_RGB, GL_UNSIGNED_BYTE, imageData);
        fbo.unbind();

        // Iterate through all the pixels in the imageBuffer and set the
        uint i,j;
        GLuint idx;
        bool flag = false;
        QImage image(w, h, QImage::Format_RGB32);
        for (i=0; i < w; ++i) {
            for (j=0; j < h; ++j) {
                uint pos = 3 * (i + j*w);
                unsigned char r = imageData[0+pos];
                unsigned char g = imageData[1+pos];
                unsigned char b = imageData[2+pos];

                // Note that cmap returns false if it is an unmapped pixel (background)
                if (cmap.colorToIdx(r,g,b,&idx)) {
                    //fprintf(stderr, "\nColor Index(%d,%d) = [%u,%u,%u] -> %u", i,j,r,g,b,idx);
                    if (idx < vvec->size()) {
                        vvec->at(idx) = true;
                    } else {
                        fprintf(stderr, "\nCritical error - vector overrun!");
                    }
                    flag = true;
                    image.setPixel(i,j,qRgb(r,g,b));
                } else {
                    image.setPixel(i,j,qRgb(0.0f,0.0f,0.0f));
                }
                //image.setPixel(i,j,qRgb(r*256.0f, g*256.0f, b*256.0f));
            }
        }
        // Write out the render zone using QImage
        image.save(QString("dump.tif"), "TIFF");

        if (!flag) fill(vvec->begin(), vvec->end(), false);
        updateGL();
    }

    // Clear away image data
    delete [] imageData;
}

