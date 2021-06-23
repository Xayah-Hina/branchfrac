#include "Common/Common.h"

#include "Demos/Common/DemoBase.h"
#include "Demos/Common/TweakBarParameters.h"
#include "Demos/Visualization/MiniGL.h"

#include "Simulation/Simulation.h"
#include "Simulation/TimeManager.h"
#include "Simulation/SimulationModel.h"
#include "Simulation/CubicSDFCollisionDetection.h"
#include "Simulation/DistanceFieldCollisionDetection.h"

#include "Utils/OBJLoader.h"
#include "Utils/SceneLoader.h"
#include "Utils/TetGenLoader.h"
#include "Utils/Logger.h"
#include "Utils/Timing.h"
#include "Utils/FileSystem.h"
#include "Utils/pymesh/MshLoader.h"

#include "GL/glut.h"

#define _USE_MATH_DEFINES

#include "math.h"

// Enable memory leak detection
#if defined(_DEBUG) && !defined(EIGEN_ALIGN)
#define new DEBUG_NEW
#endif

using namespace PBD;
using namespace Eigen;
using namespace std;
using namespace Utilities;

DemoBase *base;
Vector3r camPos;
Vector3r camLookat;
DistanceFieldCollisionDetection cd;
short simulationMethod = 2;

short clothSimulationMethod = 2;
short solidSimulationMethod = 2;
short bendingMethod = 2;
bool enableExportOBJ = false;
unsigned int exportFPS = 25;
Real nextFrameTime = 0.0;
unsigned int frameCounter = 1;

void initParameters();
void timeStep();
void buildModel();
void readScene(const bool readFile);
void render();
void reset();

void createMesh();

void seperate();


int main(int argc, char **argv)
{
    REPORT_MEMORY_LEAKS

    std::string sceneFileName;

    //    argc = 2;
    //    char *argvc[2] = {(char *)"argv1", (char *)"TreeBranchScene.json"};

    base = new DemoBase();
    base->init(argc, argv, "Test demo");

    SimulationModel *model = new SimulationModel();
    model->init();
    Simulation::getCurrent()->setModel(model);

    buildModel();

    initParameters();
    //    base->readParameters();

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
    MiniGL::setViewport(45.0f, 0.1f, 500.0f, Vector3r(0.0, 10.0, 20.0), Vector3r(0.0, 0.0, 0.0));

    glutMainLoop();

    Utilities::Timing::printAverageTimes();
    Utilities::Timing::printTimeSums();

    delete Simulation::getCurrent();
    delete base;
    delete model;

    return 0;
}

void loadbranch();
void loadObj(const std::string &filename, VertexData &vd, IndexedFaceMesh &mesh, const Vector3r &scale);
void loadfloor();

void buildModel()
{
    TimeManager::getCurrent()->setTimeStepSize(static_cast<Real>(0.005));

    SimulationModel *model = Simulation::getCurrent()->getModel();
    Simulation::getCurrent()->getTimeStep()->setCollisionDetection(*model, &cd);

//    createMesh();

    loadbranch();

    loadfloor();
//    seperate();
}

void loadbranch()
{
    vector<Vector3r> vertices;
    vector<unsigned int> tets;
    //    TetGenLoader::loadMSHModel("/Users/polaris/Desktop/tetwildsample/branch_.msh", vertices, tets);

    PyMesh::MshLoader branch("/Users/polaris/Desktop/tetwildsample/branch_.msh");

    SimulationModel *model = Simulation::getCurrent()->getModel();
    model->addTetModel((unsigned int) branch.vertices.size(), (unsigned int) branch.tets.size() / 4u, branch.vertices.data(), branch.tets.data());

    ParticleData &pd = model->getParticles();
    for (unsigned int i = 0; i < pd.getNumberOfParticles(); i++)
    {
        pd.setMass(i, 1.0);
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
            for (unsigned int i = 0; i < nTets; i++)
            {
                const unsigned int v1 = tets[4 * i];
                const unsigned int v2 = tets[4 * i + 1];
                const unsigned int v3 = tets[4 * i + 2];
                const unsigned int v4 = tets[4 * i + 3];

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
    SimulationModel::TetModelVector &tmv = model->getTetModels();
    for (unsigned int i = 0; i < tmv.size(); i++)
    {
        const unsigned int nVert = tmv[i]->getParticleMesh().numVertices();
        unsigned int offset = tmv[i]->getIndexOffset();
        tmv[i]->setFrictionCoeff(static_cast<Real>(0.1));
        cd.addCollisionObjectWithoutGeometry(i, CollisionDetection::CollisionObject::TetModelCollisionObjectType, &pd.getPosition(offset), nVert, true);
    }
}

void loadfloor()
{
    string fileName = FileSystem::normalizePath(base->getDataPath() + "/models/cube.obj");
    IndexedFaceMesh mesh;
    VertexData vd;
    loadObj(fileName, vd, mesh, Vector3r::Ones());

    SimulationModel *model = Simulation::getCurrent()->getModel();
    SimulationModel::RigidBodyVector &rb = model->getRigidBodies();
    rb.resize(1);
    rb[0] = new RigidBody();
    rb[0]->initBody(1.0, Vector3r(0.0, -5.5, 0.0), Quaternionr(1.0, 0.0, 0.0, 0.0), vd, mesh, Vector3r(100.0, 1.0, 100.0));
    rb[0]->setMass(0.0);

    Simulation::getCurrent()->getTimeStep()->setCollisionDetection(*model, &cd);
    cd.setTolerance(static_cast<Real>(0.05));

    const std::vector<Vector3r> *vertices1 = rb[0]->getGeometry().getVertexDataLocal().getVertices();
    const unsigned int nVert1 = static_cast<unsigned int>(vertices1->size());
    cd.addCollisionBox(0, CollisionDetection::CollisionObject::RigidBodyCollisionObjectType, &(*vertices1)[0], nVert1, Vector3r(100.0, 1.0, 100.0));
}

void seperate()
{
    SimulationModel *model = Simulation::getCurrent()->getModel();
    ParticleData &pd = model->getParticles();
    TetModel *tm = model->getTetModels()[0];
    TetModel::ParticleMesh &pm = tm->getParticleMesh();

    model->clearConstraint();

    int totaltet = pm.numTets();
    for (int i = 0; i < totaltet; ++i)
    {
        std::vector<unsigned int> &tetinds = pm.getTets();
        unsigned int a = tetinds[i * 4 + 0];
        unsigned int b = tetinds[i * 4 + 1];
        unsigned int c = tetinds[i * 4 + 2];
        unsigned int d = tetinds[i * 4 + 3];
        Vector3r center = pd.getPosition(a) + pd.getPosition(b) + pd.getPosition(c) + pd.getPosition(d);
        center = center / 4.0;

        unsigned int tail = pd.size();

        pd.addVertex(pd.getPosition(a));
        pd.addVertex(pd.getPosition(b));
        pd.addVertex(pd.getPosition(c));
        pd.addVertex(pd.getPosition(d));
        pd.addVertex(pd.getPosition(a));
        pd.addVertex(pd.getPosition(b));
        pd.addVertex(pd.getPosition(c));
        pd.addVertex(pd.getPosition(d));
        pd.addVertex(center);
        pd.addVertex(center);
        pd.addVertex(center);
        pd.addVertex(center);
        tetinds[i * 4 + 0] = tail + 8;

        vector<unsigned int> tet1{tail + 9, a, tail + 1, tail + 2};
        vector<unsigned int> tet2{tail + 10, tail, tail + 5, tail + 3};
        vector<unsigned int> tet3{tail + 11, tail + 4, tail + 6, tail + 7};

        pm.addTet(tet1.data());
        pm.addTet(tet2.data());
        pm.addTet(tet3.data());

        pm.addpointnum(12);
        tm->reflash();
    }

    // init constraints
    for (unsigned int cm = 0; cm < model->getTetModels().size(); cm++)
    {
        const unsigned int nTets = model->getTetModels()[cm]->getParticleMesh().numTets();
        const unsigned int *tets = model->getTetModels()[cm]->getParticleMesh().getTets().data();
        const IndexedTetMesh::VertexTets *vTets = model->getTetModels()[cm]->getParticleMesh().getVertexTets().data();

        TetModel::ParticleMesh &mesh = model->getTetModels()[cm]->getParticleMesh();
        for (unsigned int i = 0; i < nTets; i++)
        {
            const unsigned int v1 = tets[4 * i];
            const unsigned int v2 = tets[4 * i + 1];
            const unsigned int v3 = tets[4 * i + 2];
            const unsigned int v4 = tets[4 * i + 3];

            model->addFEMTetConstraint(v1, v2, v3, v4);
        }
    }
}

void createMesh()
{
    const unsigned int pointnum = 4;
    Vector3r points[pointnum];
    points[0] = Vector3r(-1, -1, 1);
    points[1] = Vector3r(1, -1, 1);
    points[2] = Vector3r(0, -1, -1);
    points[3] = Vector3r(0, 1, 1);
    //    points[4] = Vector3r(1, 1, 0);

    vector<unsigned int> indices;
    indices.push_back(0);
    indices.push_back(1);
    indices.push_back(2);
    indices.push_back(3);
    //    indices.push_back(1);
    //    indices.push_back(2);
    //    indices.push_back(3);
    //    indices.push_back(4);

    SimulationModel *model = Simulation::getCurrent()->getModel();
    model->addTetModel(pointnum, (unsigned int) indices.size() / 4u, points, indices.data());

    ParticleData &pd = model->getParticles();
    for (unsigned int i = 0; i < pd.getNumberOfParticles(); i++)
    {
        pd.setMass(i, 1.0);
    }
//    pd.setMass(0, 0.0);
    //    pd.setMass(1, 0.0);
    //    pd.setMass(2, 0.0);

    // init constraints
    for (unsigned int cm = 0; cm < model->getTetModels().size(); cm++)
    {
        const unsigned int nTets = model->getTetModels()[cm]->getParticleMesh().numTets();
        const unsigned int *tets = model->getTetModels()[cm]->getParticleMesh().getTets().data();
        const IndexedTetMesh::VertexTets *vTets = model->getTetModels()[cm]->getParticleMesh().getVertexTets().data();

        TetModel::ParticleMesh &mesh = model->getTetModels()[cm]->getParticleMesh();
        for (unsigned int i = 0; i < nTets; i++)
        {
            const unsigned int v1 = tets[4 * i];
            const unsigned int v2 = tets[4 * i + 1];
            const unsigned int v3 = tets[4 * i + 2];
            const unsigned int v4 = tets[4 * i + 3];

            model->addFEMTetConstraint(v1, v2, v3, v4);
        }
    }


    SimulationModel::TetModelVector &tmv = model->getTetModels();
    for (unsigned int i = 0; i < tmv.size(); i++)
    {
        const unsigned int nVert = tmv[i]->getParticleMesh().numVertices();
        unsigned int offset = tmv[i]->getIndexOffset();
        tmv[i]->setFrictionCoeff(static_cast<Real>(0.1));
        cd.addCollisionObjectWithoutGeometry(i, CollisionDetection::CollisionObject::TetModelCollisionObjectType, &pd.getPosition(offset), nVert, true);
    }

    LOG_INFO << "Number of tets: " << indices.size() / 4;
    LOG_INFO << "Number of vertices: " << 4;
}

// ===========================================================================================================================================
// ===========================================================================================================================================
// ===========================================================================================================================================
// ===========================================================================================================================================
// ===========================================================================================================================================
// ===========================================================================================================================================
// ===========================================================================================================================================
// ===========================================================================================================================================

void initParameters()
{
    TwRemoveAllVars(MiniGL::getTweakBar());
    TweakBarParameters::cleanup();

    MiniGL::initTweakBarParameters();

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
    }

    // Update visualization models
    const ParticleData &pd = model->getParticles();
    for (unsigned int i = 0; i < model->getTetModels().size(); i++)
    {
        model->getTetModels()[i]->updateMeshNormals(pd);
        model->getTetModels()[i]->updateVisMesh(pd);
    }
    for (unsigned int i = 0; i < model->getTriangleModels().size(); i++)
    {
        model->getTriangleModels()[i]->updateMeshNormals(pd);
    }
}

void render()
{
    base->render();
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