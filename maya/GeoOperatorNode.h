//
//  GeoOperatorNode.h
//
//  Created by Xian Xiao on 6/02/16.
//
//

#ifndef __GeoOperatorNode__
#define __GeoOperatorNode__

#include <maya/MPxLocatorNode.h>
#include <maya/MPointArray.h>
#include <maya/MBoundingBox.h>
#include <maya/MFloatPoint.h>
#include <maya/MDagPath.h>
#include <maya/MItMeshVertex.h>

#include <vector>

class GeoOperatorNode : public MPxLocatorNode
{
public:
    GeoOperatorNode();
    
    virtual ~GeoOperatorNode();
    
    virtual MStatus compute(const MPlug &plug, MDataBlock &dataBlock);
    virtual void draw(M3dView &view, const MDagPath &path, M3dView::DisplayStyle style, M3dView::DisplayStatus);
    
    static void* creator();
    static MStatus initialize();
    
    static MTypeId typeId;
    
    static MObject ia_mesh;
    static MObject ia_type;
    static MObject oa_mesh;
    
    enum { kMean, kGaussian, kMaxPrinciple, kMinPrinciple, kPricipleRatio, kErrorQuadratic, kHybrid } CurvatureType;
    
    struct Plane
    {
        Plane():a(0),b(0),c(0),d(1.){}
        Plane( const double a, const double b, const double c ) : a(a), b(b), c(c), d(1.0){}
        
        double a, b, c, d;
    };
    
private:
    
    void set_color( MDataHandle &outputHnd );
    void vertex_curvatures( MDataHandle &outputHnd, std::vector< double > &mean_curvs, short curvature_type  );
    
    double face_mix_area( MObject &mesh, MItMeshVertex &vIter ) const;
    MVector laplace_beltrami( MObject &mesh, MItMeshVertex &vIter ) const;
    
    double gaussian_curvature( MObject &mesh, MItMeshVertex &vIter ) const;
    
    bool is_vertex_obtuse( const MPoint &pt1, const MPoint &pt2, const MPoint &pt3 ) const;
    bool is_face_obtuse( const MPoint &vt, const MPoint &a, const MPoint &b ) const;
    
    void encode_color( const std::vector< double >& means, MColorArray &colors );
    
    void error_qudratic( MObject &mesh, std::vector<double> &errors );
    
    void normalize( std::vector< double > &values );
    
    
    
private:
    
    short m_type;
};

#endif
