#include "student_code.h"
#include "mutablePriorityQueue.h"

using namespace std;

namespace CGL
{

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (class member).
   *
   * @param points A vector of points in 2D
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector2D> BezierCurve::evaluateStep(std::vector<Vector2D> const &points)
  { 
    // TODO Part 1.
    vector<Vector2D> output;
    for (unsigned int i = 0; i < points.size()-1; i++) {
        Vector2D interp = points.at(i) * (1-t) + points.at(i+1) * t;
        output.push_back(interp);
    }
    return output;
  }

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (function parameter).
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector3D> BezierPatch::evaluateStep(std::vector<Vector3D> const &points, double t) const
  {
    // TODO Part 2.
      vector<Vector3D> output;
      for (unsigned int i = 0; i < points.size()-1; i++) {
          Vector3D interp = points.at(i) * (1-t) + points.at(i+1) * t;
          output.push_back(interp);
      }
      return output;
  }

  /**
   * Fully evaluates de Casteljau's algorithm for a vector of points at scalar parameter t
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate1D(std::vector<Vector3D> const &points, double t) const
  {
    // TODO Part 2.
      vector<Vector3D> interpd_pts = points;
      while (interpd_pts.size() > 1) {
          interpd_pts = evaluateStep(interpd_pts, t);
      }
      return interpd_pts.at(0);
  }

  /**
   * Evaluates the Bezier patch at parameter (u, v)
   *
   * @param u         Scalar interpolation parameter
   * @param v         Scalar interpolation parameter (along the other axis)
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate(double u, double v) const 
  {  
    // TODO Part 2.
      vector<Vector3D> interpd_pts;
      for (unsigned int i = 0; i < controlPoints.size(); i++) {
          interpd_pts.push_back(evaluate1D(controlPoints.at(i), u));
      }
      return evaluate1D(interpd_pts, v);
  }

  Vector3D Vertex::normal( void ) const
  {
    // TODO Part 3.
    // Returns an approximate unit normal at this vertex, computed by
    // taking the area-weighted average of the normals of neighboring
    // triangles, then normalizing.
      Vector3D output = Vector3D();
      HalfedgeCIter h = halfedge();
      do {
          Vector3D v1 = h->vertex()->position,
                  v2 = h->next()->vertex()->position,
                  v3 = h->next()->next()->vertex()->position;
          double area = cross(v2-v1, v3-v1).norm() / 2;

          output += h->face()->normal() * area;

          h = h->twin()->next();
      } while (h != halfedge());

      output.normalize();
      return output;
  }

  EdgeIter HalfedgeMesh::flipEdge( EdgeIter e0 )
  {
    // TODO Part 4.
    // This method should flip the given edge and return an iterator to the flipped edge.

    // Never flip a boundary edge
    if (e0->isBoundary()) {
        return e0;
    }

    // Get all the components
    HalfedgeIter bc = e0->halfedge();
    HalfedgeIter ca = bc->next();
    HalfedgeIter ab = ca->next();

    HalfedgeIter cb = bc->twin();
    HalfedgeIter bd = cb->next();
    HalfedgeIter dc = bd->next();

    FaceIter abc = bc->face();
    FaceIter bcd = cb->face();

    VertexIter a = ab->vertex();
    VertexIter b = bc->vertex();
    VertexIter c = ca->vertex();
    VertexIter d = dc->vertex();

    // Reassign the components correspondingly
    // triangle ABC -> ADC
    bc->setNeighbors(dc, cb, a, e0, abc);
    dc->setNeighbors(ca, dc->twin(), d, dc->edge(), abc);
    ca->setNeighbors(bc, ca->twin(), c, ca->edge(), abc);
    abc->halfedge() = bc;

    //triangle CBD -> DAB
    cb->setNeighbors(ab, bc, d, e0, bcd);
    ab->setNeighbors(bd, ab->twin(), a, ab->edge(), bcd);
    bd->setNeighbors(cb, bd->twin(), b, bd->edge(), bcd);
    bcd->halfedge() = cb;

    a->halfedge() = ab;
    b->halfedge() = bd;
    c->halfedge() = ca;
    d->halfedge() = dc;

    return e0;
  }

  VertexIter HalfedgeMesh::splitEdge( EdgeIter e0 )
  {
    // TODO Part 5.
    // This method should split the given edge and return an iterator to the newly inserted vertex.
    // The halfedge of this vertex should point along the edge that was split, rather than the new edges.

    // Ignore requests to split boundary edges
    if (e0->isBoundary()) {
        return e0->halfedge()->vertex();
    }

    // Get all the components
    HalfedgeIter bm = e0->halfedge();
    HalfedgeIter ca = bm->next();
    HalfedgeIter ab = ca->next();

    HalfedgeIter cm = bm->twin();
    HalfedgeIter bd = cm->next();
    HalfedgeIter dc = bd->next();

    FaceIter abm = bm->face();
    FaceIter dcm = cm->face();

    VertexIter a = ab->vertex();
    VertexIter b = bd->vertex();
    VertexIter c = ca->vertex();
    VertexIter d = dc->vertex();

    // Make new components
    VertexIter m = newVertex();
    m->position = 0.5 * (b->position + c->position);

    HalfedgeIter ma = newHalfedge(),
        am = newHalfedge(),
        mb = newHalfedge(),
        mc = newHalfedge(),
        md = newHalfedge(),
        dm = newHalfedge();

    EdgeIter amma = newEdge(),
        bmmb = newEdge(),
        cmmc = newEdge(),
        dmmd = newEdge();

    FaceIter mca = newFace(),
        mbd = newFace();

    // Reassign the components accordingly
    //triangle ABC->ACM
    ab->setNeighbors(bm, ab->twin(), a, ab->edge(), abm);
    bm->setNeighbors(ma, mb, b, bmmb, abm);
    ma->setNeighbors(ab, am, m, amma, abm);
    abm->halfedge() = ab;

    //triangle MCA
    mc->setNeighbors(ca, cm, m, cmmc, mca);
    ca->setNeighbors(am, ca->twin(), c, ca->edge(), mca);
    am->setNeighbors(mc, ma, a, amma, mca);
    mca->halfedge() = mc;

    //triangle DCB->DCM
    dc->setNeighbors(cm, dc->twin(), d, dc->edge(), dcm);
    cm->setNeighbors(md, mc, c, cmmc, dcm);
    md->setNeighbors(dc, dm, m, dmmd, dcm);
    dcm->halfedge() = dc;

    //triangle MBD
    mb->setNeighbors(bd, bm, m, bmmb, mbd);
    bd->setNeighbors(dm, bd->twin(), b, bd->edge(), mbd);
    dm->setNeighbors(mb, md, d, dmmd, mbd);
    mbd->halfedge() = mb;

    amma->halfedge() = am;
    bmmb->halfedge() = bm;
    cmmc->halfedge() = cm;
    dmmd->halfedge() = dm;

    a->halfedge() = am;
    b->halfedge() = bm;
    c->halfedge() = cm;
    d->halfedge() = dm;
    m->halfedge() = ma;

    // for loop subdivision
    amma->isNew = true;
    bmmb->isNew = false;
    cmmc->isNew = false;
    dmmd->isNew = true;
    m->isNew = true;

    return m;
  }



  void MeshResampler::upsample( HalfedgeMesh& mesh )
  {
    // TODO Part 6.
    // This routine should increase the number of triangles in the mesh using Loop subdivision.
    // One possible solution is to break up the method as listed below.

    // 1. Compute new positions for all the vertices in the input mesh, using the Loop subdivision rule,
    // and store them in Vertex::newPosition. At this point, we also want to mark each vertex as being
    // a vertex of the original mesh.

    for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
        //calculate degree and neighbor sum
        int n = 0;
        Vector3D original_neighbor_position_sum = Vector3D();
        HalfedgeCIter h = v->halfedge();
        do {
            n++;
            original_neighbor_position_sum += h->twin()->vertex()->position;
            h = h->twin()->next();
        } while (h != v->halfedge());

        double u = n == 3 ? 3./16 : 3./ (8. * n);

        v->newPosition = (1 - n * u) * v->position + u * original_neighbor_position_sum;
        v->isNew = false;
      }
    
    // 2. Compute the updated vertex positions associated with edges, and store it in Edge::newPosition.

    for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
        HalfedgeCIter h = e->halfedge();

        VertexCIter a = h->vertex();
        VertexCIter b = h->next()->vertex();
        VertexCIter d = h->next()->next()->vertex();
        VertexCIter c = h->twin()->next()->next()->vertex();

        e->newPosition = 3./8 * (a->position + b->position) + 1./8 * (c->position + d->position);
        e->isNew = false;
    }
    
    // 3. Split every edge in the mesh, in any order. For future reference, we're also going to store some
    // information about which subdivide edges come from splitting an edge in the original mesh, and which edges
    // are new, by setting the flat Edge::isNew. Note that in this loop, we only want to iterate over edges of
    // the original mesh---otherwise, we'll end up splitting edges that we just split (and the loop will never end!)

    for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
        HalfedgeIter h = e->halfedge();
        VertexCIter a = h->vertex();
        VertexCIter b = h->next()->vertex();
        if (a->isNew || b->isNew) continue;

        VertexIter m = mesh.splitEdge(e);
        m->newPosition = e->newPosition;
    }
    
    // 4. Flip any new edge that connects an old and new vertex.

    for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
        HalfedgeIter h = e->halfedge();
        VertexCIter a = h->vertex();
        VertexCIter b = h->next()->vertex();
        if (!e->isNew || a->isNew == b->isNew) continue;

        mesh.flipEdge(e);
    }

    // 5. Copy the new vertex positions into final Vertex::position.

    for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
        v->position = v->newPosition;
    }
  }
}
