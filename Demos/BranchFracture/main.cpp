#include "Common/Common.h"
#include "Demos/Visualization/MiniGL.h"
#include "Demos/Visualization/Selection.h"
#include "GL/glut.h"
#include "Simulation/TimeManager.h"
#include <Eigen/Dense>
#include "Simulation/SimulationModel.h"
#include "Simulation/TimeStepController.h"
#include <iostream>
#include "Demos/Visualization/Visualization.h"
#include "Simulation/DistanceFieldCollisionDetection.h"
#include "Utils/OBJLoader.h"
#include "Utils/Logger.h"
#include "Utils/Timing.h"
#include "Utils/FileSystem.h"
#include "Demos/Common/DemoBase.h"
#include "Demos/Common/TweakBarParameters.h"
#include "Simulation/Simulation.h"
#include "Utils/pymesh/MshLoader.h"

#include <map>

// Enable memory leak detection
#if defined(_DEBUG) && !defined(EIGEN_ALIGN)
#define new DEBUG_NEW
#endif

using namespace PBD;
using namespace Eigen;
using namespace std;
using namespace Utilities;

void initParameters();
void timeStep();
void buildModel();
void createMesh();
void render();
void reset();
void exportOBJ();
void TW_CALL setStiffness(const void *value, void *clientData);
void TW_CALL getStiffness(void *value, void *clientData);
void TW_CALL setPoissonRatio(const void *value, void *clientData);
void TW_CALL getPoissonRatio(void *value, void *clientData);
void TW_CALL setNormalizeStretch(const void *value, void *clientData);
void TW_CALL getNormalizeStretch(void *value, void *clientData);
void TW_CALL setNormalizeShear(const void *value, void *clientData);
void TW_CALL getNormalizeShear(void *value, void *clientData);
void TW_CALL setSimulationMethod(const void *value, void *clientData);
void TW_CALL getSimulationMethod(void *value, void *clientData);
void TW_CALL setContactTolerance(const void *value, void *clientData);
void TW_CALL getContactTolerance(void *value, void *clientData);
void TW_CALL setContactStiffnessRigidBody(const void *value, void *clientData);
void TW_CALL getContactStiffnessRigidBody(void *value, void *clientData);
void TW_CALL setContactStiffnessParticleRigidBody(const void *value, void *clientData);
void TW_CALL getContactStiffnessParticleRigidBody(void *value, void *clientData);

DemoBase *base;
DistanceFieldCollisionDetection cd;

const unsigned int width = 20;
const unsigned int depth = 3;
const unsigned int height = 3;
short simulationMethod = 2;

bool enableExportOBJ = false;
unsigned int exportFPS = 25;
Real nextFrameTime = 0.0;
unsigned int frameCounter = 1;

double GPLANE = 5.3;
Vector3r eyep(0.0, 30.0, 30.0);

void mcopy(unsigned int node1, unsigned int node2, ParticleData &pd)
{
    Real mass = pd.getMass(node1);
    Real invmass = pd.getInvMass(node1);
    Vector3r x0 = pd.getPosition0(node1);
    Vector3r x = pd.getPosition(node1);
    Vector3r v = pd.getVelocity(node1);
    Vector3r a = pd.getAcceleration(node1);
    Vector3r oldx = pd.getOldPosition(node1);
    Vector3r lastx = pd.getLastPosition(node1);

    pd.setMass(node2, mass);
    pd.setPosition0(node2, x0);
    pd.setPosition(node2, x);
    pd.setVelocity(node2, v);
    pd.setAcceleration(node2, a);
    pd.setOldPosition(node2, oldx);
    pd.setLastPosition(node2, lastx);
}

void loadbranch()
{
    PyMesh::MshLoader branch("/Users/polaris/Desktop/tetwildsample/branch02_.msh");
    //    PyMesh::MshLoader branch("/Users/polaris/Desktop/tetwildsample/branch04_.msh");

    for (auto &v : branch.vertices)
    {
        v *= 10.0;
        v += Vector3r(3.0, 5.0, 0.0);
    }

    SimulationModel *model = Simulation::getCurrent()->getModel();
    model->addTetModel((unsigned int) branch.vertices.size(), (unsigned int) branch.tets.size() / 4u, branch.vertices.data(), branch.tets.data());

    LOG_INFO << "Number of tets: " << branch.tets.size() / 4;
    LOG_INFO << "Number of vertices: " << branch.vertices.size();
}

void seperate()
{
    vector<unsigned int> tet1;
    vector<unsigned int> tet2;
    vector<unsigned int> fracture;

    SimulationModel *model = Simulation::getCurrent()->getModel();
    ParticleData &pd = model->getParticles();

    TetModel::ParticleMesh &pm = model->getTetModels()[0]->getParticleMesh();

    int totaltet = pm.numTets();
    for (int i = 0; i < totaltet; ++i)
    {
        std::vector<unsigned int> &tetinds = pm.getTets();
        unsigned int a = tetinds[i * 4 + 0];
        unsigned int b = tetinds[i * 4 + 1];
        unsigned int c = tetinds[i * 4 + 2];
        unsigned int d = tetinds[i * 4 + 3];
        Vector3r apos = pd.getPosition(a);
        Vector3r bpos = pd.getPosition(b);
        Vector3r cpos = pd.getPosition(c);
        Vector3r dpos = pd.getPosition(d);
        Vector3r center = (1.0 / 4.0) * (apos + bpos + cpos + dpos);

        double plane = GPLANE;
        double d1 = apos.x() - plane;
        double d2 = bpos.x() - plane;
        double d3 = cpos.x() - plane;
        double d4 = dpos.x() - plane;
        if ((d1 * d2 < 0) || (d1 * d3 < 0) || (d1 * d4 < 0) || (d2 * d3 < 0) || (d2 * d4 < 0) || (d3 * d4 < 0))
        {
            unsigned int nVert;
            unsigned int offset;
            unsigned int tail;

            tail = pd.size();
            vector<Vector3r> p1{center, pd.getPosition(b), pd.getPosition(c), pd.getPosition(d)};
            vector<unsigned int> f1{0, 1, 2, 3};
            model->addTetModel(4, (unsigned int) f1.size() / 4u, p1.data(), f1.data());
            mcopy(a, tail, pd);
            mcopy(b, tail + 1, pd);
            mcopy(c, tail + 2, pd);
            mcopy(d, tail + 3, pd);

            //            nVert = model->getTetModels().back()->getParticleMesh().numVertices();
            //            offset = model->getTetModels().back()->getIndexOffset();
            //            model->getTetModels().back()->setFrictionCoeff(static_cast<Real>(0.1));
            //            cd.addCollisionObjectWithoutGeometry(model->getTetModels().size() - 1, CollisionDetection::CollisionObject::TetModelCollisionObjectType,
            //                                                 &pd.getPosition(offset), nVert, true);

            tail = pd.size();
            vector<Vector3r> p2{pd.getPosition(a), center, pd.getPosition(c), pd.getPosition(d)};
            vector<unsigned int> f2{0, 1, 2, 3};
            model->addTetModel(4, (unsigned int) f2.size() / 4u, p2.data(), f2.data());
            mcopy(a, tail, pd);
            mcopy(b, tail + 1, pd);
            mcopy(c, tail + 2, pd);
            mcopy(d, tail + 3, pd);

            //            nVert = model->getTetModels().back()->getParticleMesh().numVertices();
            //            offset = model->getTetModels().back()->getIndexOffset();
            //            model->getTetModels().back()->setFrictionCoeff(static_cast<Real>(0.1));
            //            cd.addCollisionObjectWithoutGeometry(model->getTetModels().size() - 1, CollisionDetection::CollisionObject::TetModelCollisionObjectType,
            //                                                 &pd.getPosition(offset), nVert, true);

            tail = pd.size();
            vector<Vector3r> p3{pd.getPosition(a), pd.getPosition(b), center, pd.getPosition(d)};
            vector<unsigned int> f3{0, 1, 2, 3};
            model->addTetModel(4, (unsigned int) f3.size() / 4u, p3.data(), f3.data());
            mcopy(a, tail, pd);
            mcopy(b, tail + 1, pd);
            mcopy(c, tail + 2, pd);
            mcopy(d, tail + 3, pd);

            //            nVert = model->getTetModels().back()->getParticleMesh().numVertices();
            //            offset = model->getTetModels().back()->getIndexOffset();
            //            model->getTetModels().back()->setFrictionCoeff(static_cast<Real>(0.1));
            //            cd.addCollisionObjectWithoutGeometry(model->getTetModels().size() - 1, CollisionDetection::CollisionObject::TetModelCollisionObjectType,
            //                                                 &pd.getPosition(offset), nVert, true);

            tail = pd.size();
            vector<Vector3r> p4{pd.getPosition(a), pd.getPosition(b), pd.getPosition(c), center};
            vector<unsigned int> f4{0, 1, 2, 3};
            model->addTetModel(4, (unsigned int) f4.size() / 4u, p4.data(), f4.data());
            mcopy(a, tail, pd);
            mcopy(b, tail + 1, pd);
            mcopy(c, tail + 2, pd);
            mcopy(d, tail + 3, pd);

            //            nVert = model->getTetModels().back()->getParticleMesh().numVertices();
            //            offset = model->getTetModels().back()->getIndexOffset();
            //            model->getTetModels().back()->setFrictionCoeff(static_cast<Real>(0.1));
            //            cd.addCollisionObjectWithoutGeometry(model->getTetModels().size() - 1, CollisionDetection::CollisionObject::TetModelCollisionObjectType,
            //                                                 &pd.getPosition(offset), nVert, true);
        } else
        {
            if (d1 < 0)
            {
                tet1.push_back(a);
                tet1.push_back(b);
                tet1.push_back(c);
                tet1.push_back(d);
            } else
            {
                tet2.push_back(a);
                tet2.push_back(b);
                tet2.push_back(c);
                tet2.push_back(d);
            }
        }
    }

    model->addTetModelFromTetModelWithNewTets(Simulation::getCurrent()->getModel()->getTetModels()[0], tet2);
    //    model->addTetModelFromTetModelWithNewTets(Simulation::getCurrent()->getModel()->getTetModels()[0], fracture);

    pm.release();
    pm.changetetind(tet1);

    for (auto &tm : model->getTetModels())
    {
        tm->reflash();
    }

    model->clearConstraint();
    //    cd.cleanup();

    // init constraints
    for (unsigned int cm = 0; cm < model->getTetModels().size(); cm++)
    {
        const unsigned int nTets = model->getTetModels()[cm]->getParticleMesh().numTets();
        const unsigned int *tets = model->getTetModels()[cm]->getParticleMesh().getTets().data();
        const IndexedTetMesh::VertexTets *vTets = model->getTetModels()[cm]->getParticleMesh().getVertexTets().data();
        if (simulationMethod == 1)
        {
            const unsigned int offset = model->getTetModels()[cm]->getIndexOffset();
            const unsigned int nEdges = model->getTetModels()[cm]->getParticleMesh().numEdges();
            const IndexedTetMesh::Edge *edges = model->getTetModels()[cm]->getParticleMesh().getEdges().data();
            for (unsigned int i = 0; i < nEdges; i++)
            {
                const unsigned int v1 = edges[i].m_vert[0] + offset;
                const unsigned int v2 = edges[i].m_vert[1] + offset;

                model->addDistanceConstraint(v1, v2);
            }

            for (unsigned int i = 0; i < nTets; i++)
            {
                const unsigned int v1 = tets[4 * i];
                const unsigned int v2 = tets[4 * i + 1];
                const unsigned int v3 = tets[4 * i + 2];
                const unsigned int v4 = tets[4 * i + 3];

                model->addVolumeConstraint(v1, v2, v3, v4);
            }
        } else if (simulationMethod == 2)
        {
            TetModel::ParticleMesh &mesh = model->getTetModels()[cm]->getParticleMesh();
            const unsigned int offset = model->getTetModels()[cm]->getIndexOffset();
            for (unsigned int i = 0; i < nTets; i++)
            {
                // @author Xayah
                const unsigned int v1 = tets[4 * i] + offset;
                const unsigned int v2 = tets[4 * i + 1] + offset;
                const unsigned int v3 = tets[4 * i + 2] + offset;
                const unsigned int v4 = tets[4 * i + 3] + offset;

                model->addFEMTetConstraint(v1, v2, v3, v4);
            }
        } else if (simulationMethod == 3)
        {
            TetModel::ParticleMesh &mesh = model->getTetModels()[cm]->getParticleMesh();
            for (unsigned int i = 0; i < nTets; i++)
            {
                const unsigned int v1 = tets[4 * i];
                const unsigned int v2 = tets[4 * i + 1];
                const unsigned int v3 = tets[4 * i + 2];
                const unsigned int v4 = tets[4 * i + 3];

                model->addStrainTetConstraint(v1, v2, v3, v4);
            }
        } else if (simulationMethod == 4)
        {
            TetModel::ParticleMesh &mesh = model->getTetModels()[cm]->getParticleMesh();
            for (unsigned int i = 0; i < nTets; i++)
            {
                const unsigned int v[4] = {tets[4 * i], tets[4 * i + 1], tets[4 * i + 2], tets[4 * i + 3]};
                // Important: Divide position correction by the number of clusters
                // which contain the vertex.
                const unsigned int nc[4] = {vTets[v[0]].m_numTets, vTets[v[1]].m_numTets, vTets[v[2]].m_numTets, vTets[v[3]].m_numTets};
                model->addShapeMatchingConstraint(4, v, nc);
            }
        }
        model->getTetModels()[cm]->updateMeshNormals(pd);
    }
}

// main
int main(int argc, char **argv)
{
    REPORT_MEMORY_LEAKS

    base = new DemoBase();
    base->init(argc, argv, "Bar collision demo");

    SimulationModel *model = new SimulationModel();
    model->init();
    Simulation::getCurrent()->setModel(model);

    buildModel();

    initParameters();

    Simulation::getCurrent()->setSimulationMethodChangedCallback([&]()
                                                                 {
                                                                     reset();
                                                                     initParameters();
                                                                     base->getSceneLoader()->readParameterObject(Simulation::getCurrent()->getTimeStep());
                                                                 });

    // OpenGL
    MiniGL::setClientIdleFunc(50, timeStep);
    MiniGL::setKeyFunc(0, 'r', reset);
    MiniGL::setKeyFunc(1, 's', seperate);
    MiniGL::setClientSceneFunc(render);
    MiniGL::setViewport(40.0f, 0.1f, 500.0f, eyep, Vector3r(5.0, 0.0, 0.0));

    TwType enumType2 = TwDefineEnum("SimulationMethodType", NULL, 0);
    TwAddVarCB(MiniGL::getTweakBar(), "SimulationMethod", enumType2, setSimulationMethod, getSimulationMethod, &simulationMethod,
               " label='Simulation method' enum='0 {None}, 1 {Volume constraints}, 2 {FEM based PBD}, 3 {Strain based dynamics (no inversion handling)}, 4 {Shape matching (no inversion handling)}' group=Simulation");
    TwAddVarCB(MiniGL::getTweakBar(), "Stiffness", TW_TYPE_REAL, setStiffness, getStiffness, model,
               " label='Stiffness'  min=0.0 step=0.1 precision=4 group='Simulation' ");
    TwAddVarCB(MiniGL::getTweakBar(), "PoissonRatio", TW_TYPE_REAL, setPoissonRatio, getPoissonRatio, model,
               " label='Poisson ratio XY'  min=0.0 step=0.1 precision=4 group='Simulation' ");
    TwAddVarCB(MiniGL::getTweakBar(), "NormalizeStretch", TW_TYPE_BOOL32, setNormalizeStretch, getNormalizeStretch, model,
               " label='Normalize stretch' group='Strain based dynamics' ");
    TwAddVarCB(MiniGL::getTweakBar(), "NormalizeShear", TW_TYPE_BOOL32, setNormalizeShear, getNormalizeShear, model,
               " label='Normalize shear' group='Strain based dynamics' ");
    TwAddVarCB(MiniGL::getTweakBar(), "ContactTolerance", TW_TYPE_REAL, setContactTolerance, getContactTolerance, &cd,
               " label='Contact tolerance'  min=0.0 step=0.001 precision=3 group=Simulation ");
    TwAddVarCB(MiniGL::getTweakBar(), "ContactStiffnessRigidBody", TW_TYPE_REAL, setContactStiffnessRigidBody, getContactStiffnessRigidBody, model,
               " label='Contact stiffness RB'  min=0.0 step=0.1 precision=2 group=Simulation ");
    TwAddVarCB(MiniGL::getTweakBar(), "ContactStiffnessParticleRigidBody", TW_TYPE_REAL, setContactStiffnessParticleRigidBody, getContactStiffnessParticleRigidBody,
               model, " label='Contact stiffness Particle-RB'  min=0.0 step=0.1 precision=2 group=Simulation ");

    glutMainLoop();

    base->cleanup();

    Utilities::Timing::printAverageTimes();
    Utilities::Timing::printTimeSums();

    delete Simulation::getCurrent();
    delete base;
    delete model;

    return 0;
}

void initParameters()
{
    TwRemoveAllVars(MiniGL::getTweakBar());
    TweakBarParameters::cleanup();

    MiniGL::initTweakBarParameters();

    TwAddVarRW(MiniGL::getTweakBar(), "ExportOBJ", TW_TYPE_BOOL32, &enableExportOBJ, " label='Export OBJ'");
    TwAddVarRW(MiniGL::getTweakBar(), "ExportFPS", TW_TYPE_UINT32, &exportFPS, " label='Export FPS'");

    TweakBarParameters::createParameterGUI();
    TweakBarParameters::createParameterObjectGUI(base);
    TweakBarParameters::createParameterObjectGUI(Simulation::getCurrent());
    TweakBarParameters::createParameterObjectGUI(Simulation::getCurrent()->getModel());
    TweakBarParameters::createParameterObjectGUI(Simulation::getCurrent()->getTimeStep());
}

void reset()
{
    Utilities::Timing::printAverageTimes();
    Utilities::Timing::reset();

    Simulation::getCurrent()->reset();
    base->getSelectedParticles().clear();

    Simulation::getCurrent()->getModel()->cleanup();
    Simulation::getCurrent()->getTimeStep()->getCollisionDetection()->cleanup();

//    nextFrameTime = 0.0;
//    frameCounter = 1;

    buildModel();
}

void timeStep()
{
    const Real pauseAt = base->getValue<Real>(DemoBase::PAUSE_AT);
    if ((pauseAt > 0.0) && (pauseAt < TimeManager::getCurrent()->getTime()))
        base->setValue(DemoBase::PAUSE, true);

    if (base->getValue<bool>(DemoBase::PAUSE))
        return;

    // Simulation code
    SimulationModel *model = Simulation::getCurrent()->getModel();
    const unsigned int numSteps = base->getValue<unsigned int>(DemoBase::NUM_STEPS_PER_RENDER);
    for (unsigned int i = 0; i < numSteps; i++)
    {
        START_TIMING("SimStep");
        Simulation::getCurrent()->getTimeStep()->step(*model);
        STOP_TIMING_AVG;

//        exportOBJ();
    }

    for (unsigned int i = 0; i < model->getTetModels().size(); i++)
        model->getTetModels()[i]->updateMeshNormals(model->getParticles());
}

void loadObj(const std::string &filename, VertexData &vd, IndexedFaceMesh &mesh, const Vector3r &scale)
{
    std::vector<OBJLoader::Vec3f> x;
    std::vector<OBJLoader::Vec3f> normals;
    std::vector<OBJLoader::Vec2f> texCoords;
    std::vector<MeshFaceIndices> faces;
    OBJLoader::Vec3f s = {(float) scale[0], (float) scale[1], (float) scale[2]};
    OBJLoader::loadObj(filename, &x, &faces, &normals, &texCoords, s);

    mesh.release();
    const unsigned int nPoints = (unsigned int) x.size();
    const unsigned int nFaces = (unsigned int) faces.size();
    const unsigned int nTexCoords = (unsigned int) texCoords.size();
    mesh.initMesh(nPoints, nFaces * 2, nFaces);
    vd.reserve(nPoints);
    for (unsigned int i = 0; i < nPoints; i++)
    {
        vd.addVertex(Vector3r(x[i][0], x[i][1], x[i][2]));
    }
    for (unsigned int i = 0; i < nTexCoords; i++)
    {
        mesh.addUV(texCoords[i][0], texCoords[i][1]);
    }
    for (unsigned int i = 0; i < nFaces; i++)
    {
        // Reduce the indices by one
        int posIndices[3];
        int texIndices[3];
        for (int j = 0; j < 3; j++)
        {
            posIndices[j] = faces[i].posIndices[j] - 1;
            if (nTexCoords > 0)
            {
                texIndices[j] = faces[i].texIndices[j] - 1;
                mesh.addUVIndex(texIndices[j]);
            }
        }

        mesh.addFace(&posIndices[0]);
    }
    mesh.buildNeighbors();

    mesh.updateNormals(vd, 0);
    mesh.updateVertexNormals(vd);

    LOG_INFO << "Number of triangles: " << nFaces;
    LOG_INFO << "Number of vertices: " << nPoints;
}

void buildModel()
{
    TimeManager::getCurrent()->setTimeStepSize(static_cast<Real>(0.005));
    SimulationModel *model = Simulation::getCurrent()->getModel();

    createMesh();

    // create static rigid body
    string fileName = FileSystem::normalizePath(base->getDataPath() + "/models/cube.obj");
    IndexedFaceMesh mesh;
    VertexData vd;
    loadObj(fileName, vd, mesh, Vector3r::Ones());

    string fileNameTorus = FileSystem::normalizePath(base->getDataPath() + "/models/torus.obj");
    IndexedFaceMesh meshTorus;
    VertexData vdTorus;
    loadObj(fileNameTorus, vdTorus, meshTorus, Vector3r::Ones());

    SimulationModel::RigidBodyVector &rb = model->getRigidBodies();
    rb.resize(1);

    // floor
    rb[0] = new RigidBody();
    rb[0]->initBody(1.0, Vector3r(0.0, -5.5, 0.0), Quaternionr(1.0, 0.0, 0.0, 0.0), vd, mesh, Vector3r(100.0, 1.0, 100.0));
    rb[0]->setMass(0.0);

//    // torus
//    rb[1] = new RigidBody();
//    rb[1]->initBody(1.0, Vector3r(5.0, -1.5, 0.0), Quaternionr(1.0, 0.0, 0.0, 0.0), vdTorus, meshTorus, Vector3r(3.0, 3.0, 3.0));
//    rb[1]->setMass(0.0);
//    rb[1]->setFrictionCoeff(static_cast<Real>(0.1));

    Simulation::getCurrent()->getTimeStep()->setCollisionDetection(*model, &cd);
    cd.setTolerance(static_cast<Real>(0.05));

    const std::vector<Vector3r> *vertices1 = rb[0]->getGeometry().getVertexDataLocal().getVertices();
    const unsigned int nVert1 = static_cast<unsigned int>(vertices1->size());
    cd.addCollisionBox(0, CollisionDetection::CollisionObject::RigidBodyCollisionObjectType, &(*vertices1)[0], nVert1, Vector3r(100.0, 1.0, 100.0));

//    const std::vector<Vector3r> *vertices2 = rb[1]->getGeometry().getVertexDataLocal().getVertices();
//    const unsigned int nVert2 = static_cast<unsigned int>(vertices2->size());
//    cd.addCollisionTorus(1, CollisionDetection::CollisionObject::RigidBodyCollisionObjectType, &(*vertices2)[0], nVert2, Vector2r(3.0, 1.5));

    SimulationModel::TetModelVector &tm = model->getTetModels();
    ParticleData &pd = model->getParticles();
    for (unsigned int i = 0; i < tm.size(); i++)
    {
        const unsigned int nVert = tm[i]->getParticleMesh().numVertices();
        unsigned int offset = tm[i]->getIndexOffset();
        tm[i]->setFrictionCoeff(static_cast<Real>(0.1));
        cd.addCollisionObjectWithoutGeometry(i, CollisionDetection::CollisionObject::TetModelCollisionObjectType, &pd.getPosition(offset), nVert, true);
    }
    //    seperate();
}

void render()
{
    base->render();
}

void createMesh()
{
    Vector3r points[width * height * depth];
    for (unsigned int i = 0; i < width; i++)
    {
        for (unsigned int j = 0; j < height; j++)
        {
            for (unsigned int k = 0; k < depth; k++)
            {
                points[i * height * depth + j * depth + k] = 0.3 * Vector3r((Real) i, (Real) j + static_cast<Real>(5.0), (Real) k);
            }
        }
    }

    vector<unsigned int> indices;
    for (unsigned int i = 0; i < width - 1; i++)
    {
        for (unsigned int j = 0; j < height - 1; j++)
        {
            for (unsigned int k = 0; k < depth - 1; k++)
            {
                // For each block, the 8 corners are numerated as:
                //     4*-----*7
                //     /|    /|
                //    / |   / |
                //  5*-----*6 |
                //   | 0*--|--*3
                //   | /   | /
                //   |/    |/
                //  1*-----*2
                unsigned int p0 = i * height * depth + j * depth + k;
                unsigned int p1 = p0 + 1;
                unsigned int p3 = (i + 1) * height * depth + j * depth + k;
                unsigned int p2 = p3 + 1;
                unsigned int p7 = (i + 1) * height * depth + (j + 1) * depth + k;
                unsigned int p6 = p7 + 1;
                unsigned int p4 = i * height * depth + (j + 1) * depth + k;
                unsigned int p5 = p4 + 1;

                // Ensure that neighboring tetras are sharing faces
                if ((i + j + k) % 2 == 1)
                {
                    indices.push_back(p2);
                    indices.push_back(p1);
                    indices.push_back(p6);
                    indices.push_back(p3);
                    indices.push_back(p6);
                    indices.push_back(p3);
                    indices.push_back(p4);
                    indices.push_back(p7);
                    indices.push_back(p4);
                    indices.push_back(p1);
                    indices.push_back(p6);
                    indices.push_back(p5);
                    indices.push_back(p3);
                    indices.push_back(p1);
                    indices.push_back(p4);
                    indices.push_back(p0);
                    indices.push_back(p6);
                    indices.push_back(p1);
                    indices.push_back(p4);
                    indices.push_back(p3);
                } else
                {
                    indices.push_back(p0);
                    indices.push_back(p2);
                    indices.push_back(p5);
                    indices.push_back(p1);
                    indices.push_back(p7);
                    indices.push_back(p2);
                    indices.push_back(p0);
                    indices.push_back(p3);
                    indices.push_back(p5);
                    indices.push_back(p2);
                    indices.push_back(p7);
                    indices.push_back(p6);
                    indices.push_back(p7);
                    indices.push_back(p0);
                    indices.push_back(p5);
                    indices.push_back(p4);
                    indices.push_back(p0);
                    indices.push_back(p2);
                    indices.push_back(p7);
                    indices.push_back(p5);
                }
            }
        }
    }
    SimulationModel *model = Simulation::getCurrent()->getModel();
    //    model->addTetModel(width * height * depth, (unsigned int) indices.size() / 4u, points, indices.data());
    //
    //    vector<Vector3r> p1{Vector3r(5, 5, 5), Vector3r(5, 5, 4), Vector3r(5, 4, 5), Vector3r(4, 5, 5)};
    //    vector<unsigned int> f1{0, 1, 2, 3};
    //    model->addTetModel(4, (unsigned int) f1.size() / 4u, p1.data(), f1.data());
    //
    //    vector<Vector3r> p2{Vector3r(7, 7, 7), Vector3r(7, 7, 6), Vector3r(7, 6, 7), Vector3r(6, 7, 7)};
    //    vector<unsigned int> f2{0, 1, 2, 3};
    //    model->addTetModel(4, (unsigned int) f2.size() / 4u, p2.data(), f2.data());

    loadbranch();

    ParticleData &pd = model->getParticles();
    for (unsigned int i = 0; i < pd.getNumberOfParticles(); i++)
    {
        pd.setMass(i, 1.0);
        if (pd.getPosition0(i).x() < GPLANE - 0.3)
            pd.setMass(i, 0.0);
    }


    // init constraints
    for (unsigned int cm = 0; cm < model->getTetModels().size(); cm++)
    {
        const unsigned int nTets = model->getTetModels()[cm]->getParticleMesh().numTets();
        const unsigned int *tets = model->getTetModels()[cm]->getParticleMesh().getTets().data();
        const IndexedTetMesh::VertexTets *vTets = model->getTetModels()[cm]->getParticleMesh().getVertexTets().data();
        if (simulationMethod == 1)
        {
            const unsigned int offset = model->getTetModels()[cm]->getIndexOffset();
            const unsigned int nEdges = model->getTetModels()[cm]->getParticleMesh().numEdges();
            const IndexedTetMesh::Edge *edges = model->getTetModels()[cm]->getParticleMesh().getEdges().data();
            for (unsigned int i = 0; i < nEdges; i++)
            {
                const unsigned int v1 = edges[i].m_vert[0] + offset;
                const unsigned int v2 = edges[i].m_vert[1] + offset;

                model->addDistanceConstraint(v1, v2);
            }

            for (unsigned int i = 0; i < nTets; i++)
            {
                const unsigned int v1 = tets[4 * i];
                const unsigned int v2 = tets[4 * i + 1];
                const unsigned int v3 = tets[4 * i + 2];
                const unsigned int v4 = tets[4 * i + 3];

                model->addVolumeConstraint(v1, v2, v3, v4);
            }
        } else if (simulationMethod == 2)
        {
            TetModel::ParticleMesh &mesh = model->getTetModels()[cm]->getParticleMesh();
            const unsigned int offset = model->getTetModels()[cm]->getIndexOffset();
            for (unsigned int i = 0; i < nTets; i++)
            {
                // @author Xayah
                const unsigned int v1 = tets[4 * i] + offset;
                const unsigned int v2 = tets[4 * i + 1] + offset;
                const unsigned int v3 = tets[4 * i + 2] + offset;
                const unsigned int v4 = tets[4 * i + 3] + offset;

                model->addFEMTetConstraint(v1, v2, v3, v4);
            }
        } else if (simulationMethod == 3)
        {
            TetModel::ParticleMesh &mesh = model->getTetModels()[cm]->getParticleMesh();
            for (unsigned int i = 0; i < nTets; i++)
            {
                const unsigned int v1 = tets[4 * i];
                const unsigned int v2 = tets[4 * i + 1];
                const unsigned int v3 = tets[4 * i + 2];
                const unsigned int v4 = tets[4 * i + 3];

                model->addStrainTetConstraint(v1, v2, v3, v4);
            }
        } else if (simulationMethod == 4)
        {
            TetModel::ParticleMesh &mesh = model->getTetModels()[cm]->getParticleMesh();
            for (unsigned int i = 0; i < nTets; i++)
            {
                const unsigned int v[4] = {tets[4 * i], tets[4 * i + 1], tets[4 * i + 2], tets[4 * i + 3]};
                // Important: Divide position correction by the number of clusters
                // which contain the vertex.
                const unsigned int nc[4] = {vTets[v[0]].m_numTets, vTets[v[1]].m_numTets, vTets[v[2]].m_numTets, vTets[v[3]].m_numTets};
                model->addShapeMatchingConstraint(4, v, nc);
            }
        }
        model->getTetModels()[cm]->updateMeshNormals(pd);
    }

    LOG_INFO << "Number of tets: " << indices.size() / 4;
    LOG_INFO << "Number of vertices: " << width * height * depth;

}

void exportMeshOBJ(const std::string &exportFileName, const unsigned int nVert, const Vector3r *pos, const unsigned int nTri, const unsigned int *faces)
{
    // Open the file
    std::ofstream outfile(exportFileName);
    if (!outfile)
    {
        LOG_WARN << "Cannot open a file to save OBJ mesh.";
        return;
    }

    // Header
    outfile << "# Created by the PositionBasedDynamics library\n";
    outfile << "g default\n";

    // Vertices
    {
        for (auto j = 0u; j < nVert; j++)
        {
            const Vector3r &x = pos[j];
            outfile << "v " << x[0] << " " << x[1] << " " << x[2] << "\n";
        }
    }

    // faces
    {
        for (auto j = 0; j < nTri; j++)
        {
            outfile << "f " << faces[3 * j + 0] + 1 << " " << faces[3 * j + 1] + 1 << " " << faces[3 * j + 2] + 1 << "\n";
        }
    }
    outfile.close();
}

void exportOBJ()
{
    if (!enableExportOBJ)
        return;

    cout << "RENDER A FRAME" << endl;

    if (TimeManager::getCurrent()->getTime() < nextFrameTime)
        return;

    nextFrameTime += 1.0 / (Real) exportFPS;

    //////////////////////////////////////////////////////////////////////////
    // rigid bodies
    //////////////////////////////////////////////////////////////////////////

    //    std::string exportPath = base->getOutputPath() + "/export";
    //    FileSystem::makeDirs(exportPath);

    std::string exportPath = "/Users/polaris/Desktop/objs/branchfracture";
    FileSystem::makeDirs(exportPath);

    SimulationModel *model = Simulation::getCurrent()->getModel();
    const ParticleData &pd = model->getParticles();
    for (unsigned int i = 0; i < model->getTriangleModels().size(); i++)
    {
        const IndexedFaceMesh &mesh = model->getTriangleModels()[i]->getParticleMesh();
        const unsigned int offset = model->getTriangleModels()[i]->getIndexOffset();
        const Vector3r *x = model->getParticles().getVertices()->data();

        std::string fileName = "triangle_model";
        fileName = fileName + std::to_string(i) + "_" + std::to_string(frameCounter) + ".obj";
        std::string exportFileName = FileSystem::normalizePath(exportPath + "/" + fileName);

        exportMeshOBJ(exportFileName, mesh.numVertices(), &x[offset], mesh.numFaces(), mesh.getFaces().data());
    }

    for (unsigned int i = 0; i < model->getTetModels().size(); i++)
    {
        // @author Xayah
        const IndexedFaceMesh &mesh = model->getTetModels()[i]->getSurfaceMesh();
        const unsigned int offset = model->getTetModels()[i]->getIndexOffset();
        const Vector3r *x = model->getTetModels()[i]->getVisVertices().getVertices()->data();

        std::vector<unsigned int> faces = model->getTetModels()[i]->getSurfaceMesh().getFaces();
        map<unsigned int, unsigned int> mymap;

        for (int j = 0; j < faces.size(); ++j)
        {
            ++mymap[faces[j]];
        }

        vector<Vector3r> vertices;
        std::vector<unsigned int> indices;

        int num = 1;
        for (auto it = mymap.begin(); it != mymap.end(); ++it)
        {
            vertices.push_back(pd.getPosition(offset + (it->first)));
            it->second = num++;
        }

        for (int j = 0; j < faces.size(); ++j)
        {
            indices.push_back(mymap[faces[j]]);
        }

        std::string fileName = "tet_model";
        fileName = fileName + std::to_string(i) + "_" + std::to_string(frameCounter) + ".obj";
        std::string exportFileName = FileSystem::normalizePath(exportPath + "/" + fileName);

        exportMeshOBJ(exportFileName, vertices.size(), vertices.data(), indices.size() / 3u, indices.data());
    }

    for (unsigned int i = 0; i < model->getRigidBodies().size(); i++)
    {
        const IndexedFaceMesh &mesh = model->getRigidBodies()[i]->getGeometry().getMesh();
        const Vector3r *x = model->getRigidBodies()[i]->getGeometry().getVertexData().getVertices()->data();

        std::string fileName = "rigid_body";
        fileName = fileName + std::to_string(i) + "_" + std::to_string(frameCounter) + ".obj";
        std::string exportFileName = FileSystem::normalizePath(exportPath + "/" + fileName);

        exportMeshOBJ(exportFileName, mesh.numVertices(), x, mesh.numFaces(), mesh.getFaces().data());
    }

    frameCounter++;
}

void TW_CALL setStiffness(const void *value, void *clientData)
{
    const Real val = *(const Real *) (value);
    ((SimulationModel *) clientData)->setValue<Real>(SimulationModel::SOLID_STIFFNESS, val);
}

void TW_CALL getStiffness(void *value, void *clientData)
{
    *(Real *) (value) = ((SimulationModel *) clientData)->getValue<Real>(SimulationModel::SOLID_STIFFNESS);
}

void TW_CALL setPoissonRatio(const void *value, void *clientData)
{
    const Real val = *(const Real *) (value);
    ((SimulationModel *) clientData)->setValue<Real>(SimulationModel::SOLID_POISSON_RATIO, val);
}

void TW_CALL getPoissonRatio(void *value, void *clientData)
{
    *(Real *) (value) = ((SimulationModel *) clientData)->getValue<Real>(SimulationModel::SOLID_POISSON_RATIO);
}

void TW_CALL setNormalizeStretch(const void *value, void *clientData)
{
    const bool val = *(const bool *) (value);
    ((SimulationModel *) clientData)->setValue<Real>(SimulationModel::SOLID_NORMALIZE_STRETCH, val);
}

void TW_CALL getNormalizeStretch(void *value, void *clientData)
{
    *(bool *) (value) = ((SimulationModel *) clientData)->getValue<Real>(SimulationModel::SOLID_NORMALIZE_STRETCH);
}

void TW_CALL setNormalizeShear(const void *value, void *clientData)
{
    const bool val = *(const bool *) (value);
    ((SimulationModel *) clientData)->setValue<Real>(SimulationModel::SOLID_NORMALIZE_SHEAR, val);
}

void TW_CALL getNormalizeShear(void *value, void *clientData)
{
    *(bool *) (value) = ((SimulationModel *) clientData)->getValue<Real>(SimulationModel::SOLID_NORMALIZE_SHEAR);
}

void TW_CALL setSimulationMethod(const void *value, void *clientData)
{
    const short val = *(const short *) (value);
    *((short *) clientData) = val;
    reset();
}

void TW_CALL getSimulationMethod(void *value, void *clientData)
{
    *(short *) (value) = *((short *) clientData);
}

void TW_CALL setContactStiffnessRigidBody(const void *value, void *clientData)
{
    const Real val = *(const Real *) (value);
    ((SimulationModel *) clientData)->setContactStiffnessRigidBody(val);
}

void TW_CALL getContactStiffnessRigidBody(void *value, void *clientData)
{
    *(Real *) (value) = ((SimulationModel *) clientData)->getContactStiffnessRigidBody();
}

void TW_CALL setContactStiffnessParticleRigidBody(const void *value, void *clientData)
{
    const Real val = *(const Real *) (value);
    ((SimulationModel *) clientData)->setContactStiffnessParticleRigidBody(val);
}

void TW_CALL getContactStiffnessParticleRigidBody(void *value, void *clientData)
{
    *(Real *) (value) = ((SimulationModel *) clientData)->getContactStiffnessParticleRigidBody();
}

void TW_CALL setContactTolerance(const void *value, void *clientData)
{
    const Real val = *(const Real *) (value);
    ((DistanceFieldCollisionDetection *) clientData)->setTolerance(val);
}

void TW_CALL getContactTolerance(void *value, void *clientData)
{
    *(Real *) (value) = ((DistanceFieldCollisionDetection *) clientData)->getTolerance();
}
