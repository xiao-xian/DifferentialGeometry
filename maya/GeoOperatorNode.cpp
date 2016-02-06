//
//  GeoOperatorNode.cpp
//
//  Created by Xian Xiao on 6/02/16.
//
//

#include "GeoOperatorNode.h"

#include <maya/MFnCompoundAttribute.h>
#include <maya/MFnStringData.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnPluginData.h>
#include <maya/MFnEnumAttribute.h>
#include <maya/MArrayDataHandle.h>
#include <maya/MPlugArray.h>
#include <maya/MDagPath.h>
#include <maya/MFnEnumAttribute.h>
#include <maya/MFloatPointArray.h>
#include <maya/MMeshIntersector.h>
#include <maya/MFnCamera.h>
#include <maya/MFnMesh.h>
#include <maya/MTime.h>
#include <maya/MAnimControl.h>
#include <maya/MFnMeshData.h>
#include <maya/MItMeshEdge.h>
#include <maya/MItMeshVertex.h>
#include <maya/MItMeshFaceVertex.h>
#include <maya/MItMeshPolygon.h>

#include <maya/MMeshIntersector.h>


MObject GeoOperatorNode::ia_mesh;
MObject GeoOperatorNode::ia_type;
MObject GeoOperatorNode::oa_mesh;


MTypeId GeoOperatorNode::typeId(0x001151AC);

GeoOperatorNode::GeoOperatorNode() : m_type( 0 )
{
}

GeoOperatorNode::~GeoOperatorNode()
{
}

void* GeoOperatorNode::creator()
{
    return new GeoOperatorNode;
}

MStatus GeoOperatorNode::initialize()
{
    MStatus stat;
    
    {
        MFnTypedAttribute tAttr;
        ia_mesh = tAttr.create( "inputMesh", "inputMesh", MFnData::kMesh, MObject::kNullObj, &stat );
        tAttr.setDisconnectBehavior( MFnAttribute::kDelete );
        CHECK_MSTATUS( stat );
        addAttribute( ia_mesh );
    }
    
    {
        MFnEnumAttribute eAttr;
        ia_type = eAttr.create( "Type", "Type" );
        eAttr.addField( "Mean", 0 );
        eAttr.addField( "Gaussian" , 1);
        eAttr.addField( "MaxPrinciple", 2 );
        eAttr.addField( "MinPrinciple", 3 );
        addAttribute( ia_type );
    }
    
    {
        MFnTypedAttribute tAttr;
        oa_mesh = tAttr.create( "outputMesh", "outputMesh", MFnData::kMesh, MObject::kNullObj, &stat );
        tAttr.setDisconnectBehavior( MFnAttribute::kDelete );
        CHECK_MSTATUS( stat );
        addAttribute( oa_mesh );
    }
    
    attributeAffects( ia_mesh, oa_mesh );
    attributeAffects( ia_type, oa_mesh );
    return stat;
}

MStatus GeoOperatorNode::compute(const MPlug &plug, MDataBlock &dataBlock)
{
    if( plug == oa_mesh )
    {
        MDataHandle inputHand = dataBlock.inputValue( ia_mesh );
        m_type = dataBlock.inputValue( ia_type ).asShort();
        MDataHandle outputHnd = dataBlock.outputValue( oa_mesh );
        outputHnd.copy( inputHand );
        
        set_color( outputHnd );
        
        return dataBlock.setClean( plug );
    }
    return MS::kUnknownParameter;
}

void GeoOperatorNode::draw(M3dView &view, const MDagPath &path, M3dView::DisplayStyle style, M3dView::DisplayStatus)
{
    
}

void GeoOperatorNode::set_color( MDataHandle &outputHnd )
{
    std::vector< double > means;
    vertex_curvatures( outputHnd, means );
    MStatus stat;
    MFnMesh meshFn( outputHnd.asMesh(), &stat );
    
    MString color_set("MyColor");
    meshFn.createColorSetDataMesh( color_set );
    MColorArray vertex_colors;
    MPointArray pts;
    meshFn.getPoints( pts );
    
    MColorArray curvature_colors;
    encode_color( means,  curvature_colors );
    
    MIntArray poly_vertId;
    for( size_t f = 0; f < meshFn.numPolygons(); f++ )
    {
        MIntArray vertexList;
        meshFn.getPolygonVertices( f, vertexList );
        for( size_t v = 0; v < vertexList.length(); v++ )
        {
            poly_vertId.append( vertex_colors.length() );
            vertex_colors.append( curvature_colors[ vertexList[ v ] ] );
        }
    }
    
    stat = meshFn.setColors( vertex_colors, &color_set );
    stat = meshFn.assignColors( poly_vertId, &color_set );
}

void GeoOperatorNode::vertex_curvatures( MDataHandle &outputHnd, std::vector< double > &curvatures )
{
    MObject meshData = outputHnd.asMesh();
    MItMeshVertex vertIter( meshData );
    MItMeshPolygon polyIter( meshData );
    
    MFnMesh meshFn( meshData );
    MPointArray pts;
    meshFn.getPoints( pts );
    curvatures.resize( pts.length() );
    
    int id = 0;
    for( ; !vertIter.isDone(); vertIter.next() )
    {
        const double area = face_mix_area( meshData, vertIter );
        MVector normal = area > 0.0 ? ( 1./area ) * laplace_beltrami( meshData, vertIter ) : MVector(0, 0, 0 );
        
        const double Kh = 0.5 * normal.length();
        const double Kg = gaussian_curvature( meshData,  vertIter );
        const double deltaX = Kh*Kh - Kg;
        
        if( m_type == kMean )
        {
            curvatures[ id++ ] = Kh;
        }
        else if( m_type == kGaussian )
        {
            curvatures[ id++ ] = Kg;
        }
        else if( m_type == kMaxPrinciple )
        {
            curvatures[ id++ ] = Kh + sqrt( deltaX );
        }
        else if( m_type == kMinPrinciple )
        {
            curvatures[ id++ ] = Kh - sqrt( deltaX );
        }
    }
    
    double max_c = 0.0, min_c = 1.e+12;
    
    for( size_t i = 0; i < curvatures.size(); i++ )
    {
        max_c = ( curvatures[ i ] > max_c ? curvatures[ i ] : max_c );
        min_c = ( curvatures[ i ] < min_c ? curvatures[ i ] : min_c );
    }
    
    if( max_c > min_c )
    {
        for( size_t i = 0; i < curvatures.size(); i++ )
        {
            curvatures[ i ] = ( curvatures[ i ] - min_c ) / ( max_c - min_c );
        }
    }
}

double GeoOperatorNode::face_mix_area( MObject &mesh, MItMeshVertex &vIter ) const
{
    const MPoint current = vIter.position();
    MItMeshPolygon polyIter( mesh );
    
    MIntArray faceList;
    vIter.getConnectedFaces( faceList );
    
    int preId;
    double mix_area = 0;
    for( size_t f = 0; f < faceList.length(); f++ )
    {
        polyIter.setIndex( faceList[ f ], preId );
        
        MPointArray vertices;
        polyIter.getPoints( vertices );
        
        assert( vertices.length() == 3 );
        
        double area;
        polyIter.getArea( area );
       
        MPointArray opposites;
        for( size_t vv = 0; vv < 3; vv++ )
        {
            if( vertices[ vv ].distanceTo( current ) > 1.e-3 )
            {
                opposites.append( vertices[ vv ] );
            }
        }
        
        if( is_face_obtuse( vertices[0], vertices[1], vertices[2] ) )
        {
            if( is_vertex_obtuse( current, opposites[0], opposites[1] ) )
            {
                mix_area += 0.5 * area;
            }
            else
            {
                mix_area += 0.25 * area;
            }
        }
        else
        {
            MVector ab = current - opposites[0];
            MVector cb = opposites[1] - opposites[0];
            
            MVector ac = current - opposites[1];
            MVector bc = opposites[ 0 ] - opposites[1];
            

            ab.normalize();
            cb.normalize();
            ac.normalize();
            bc.normalize();
            
            const double angle_a = acos( ab * cb ), angle_b = acos( ac * bc );
            const double cot_a = ( angle_a != 0 ? 1./tan( angle_a ) : 0 );
            const double cot_b = ( angle_b != 0 ? 1./tan( angle_b ) : 0 );
            
            const double la = current.distanceTo( opposites[ 0 ] );
            const double lb = current.distanceTo( opposites[ 1 ] );
            
            mix_area += 0.125 * ( cot_b * la*la + cot_a *lb * lb );
        }
    }
    
    return mix_area;
    
}

double GeoOperatorNode::gaussian_curvature( MObject &mesh, MItMeshVertex &vIter ) const
{
    const MPoint current = vIter.position();
    MItMeshPolygon polyIter( mesh );
    
    MIntArray faceList;
    vIter.getConnectedFaces( faceList );
    
    int preId;
    double sum_theta = 0.;
    for( size_t f = 0; f < faceList.length(); f++ )
    {
        polyIter.setIndex( faceList[ f ], preId );
        
        MPointArray vertices;
        polyIter.getPoints( vertices );
        
        assert( vertices.length() == 3 );
        
        MPointArray opposites;
        for( size_t vv = 0; vv < 3; vv++ )
        {
            if( vertices[ vv ].distanceTo( current ) > 1.e-3 )
            {
                opposites.append( vertices[ vv ] );
            }
        }
        
        MVector ba = opposites[ 0 ] - current;
        MVector ca = opposites[ 1 ] - current;
        ba.normalize();
        ca.normalize();
        
        sum_theta += acos( ba * ca );
    }
    
    double mix_area = face_mix_area( mesh, vIter );
    
    return ( 2.* M_PI - sum_theta ) / mix_area;
}

bool GeoOperatorNode::is_face_obtuse( const MPoint &p1, const MPoint &p2, const MPoint &p3 ) const
{
    return ( is_vertex_obtuse( p1, p2, p3 ) || is_vertex_obtuse( p2, p1, p3 ) || is_vertex_obtuse( p3, p1, p2 ) );
}

bool GeoOperatorNode::is_vertex_obtuse( const MPoint &vt, const MPoint &a, const MPoint &b ) const
{
    const double ab = (a - b ).length();
    const double va = (vt - a ).length();
    const double vb = (vt - b ).length();
    
    return va*va + vb*vb < ab*ab;
}

MVector GeoOperatorNode::laplace_beltrami( MObject &mesh, MItMeshVertex &vIter )  const
{
    const MPoint current = vIter.position();
    
    MItMeshEdge edgeIter( mesh );
    MItMeshPolygon polyIter( mesh );
   
    MIntArray edgeList;
    vIter.getConnectedEdges( edgeList );
    
    int preId;
    MVector normal( 0,0,0);
    for( unsigned int e = 0; e < edgeList.length(); e++ )
    {
        edgeIter.setIndex( edgeList[e], preId );
        
        MPoint pt1 = edgeIter.point( 0 ), pt2 = edgeIter.point( 1 );
        
        MIntArray faceList;
        edgeIter.getConnectedFaces( faceList );
        
        int preFaceId;
        double cot_w = 0.;

        {
            for( size_t f = 0; f < faceList.length(); f++ )
            {
                polyIter.setIndex( faceList[ f ], preFaceId );
                
                MPointArray vertices;
                polyIter.getPoints( vertices );
                
                assert( vertices.length() == 3 );
                
                for( size_t v = 0; v < 3; v++ )
                {
                    if( vertices[ v ].distanceTo( pt1 ) > 1.e-3 && vertices[ v ].distanceTo( pt2 ) > 1.e-3 )
                    {
                        MVector ba = pt1 - vertices[ v ];
                        MVector ca = pt2 - vertices[ v ];
                        
                        ba.normalize();
                        ca.normalize();
                        
                        double angle = acos( ba*ca );
                        cot_w += ( angle != 0 ? 1./tan( angle ) : 0.0 );
                    }
                }
            }
        }
        
        MPoint the_other = ( current.distanceTo( pt1 ) < 1.e-3 ? pt2 : pt1 );
        
        normal += cot_w * ( current - the_other );
    }

    return normal;
}

void GeoOperatorNode::encode_color( const std::vector< double >& means, MColorArray &colors )
{
    const MColor red( 1.0, 0.0, 0.0 ), yellow( 1.0, 1.0, 0.0 ), green( 0.0, 1.0, 0.0 ), cyan( 0, 1., 1 ), blue( 0.0, 0.0, 1.0 );
    const double inter[5] = { 1.0, 0.75, 0.5, 0.25, 0.0 };
    
    const double sigma = 0.1;
    colors.setLength( means.size() );
    for( size_t i = 0; i < means.size(); i++ )
    {
        const double value = means[ i ];
        
        double weights[5];
        double sum = 0;
        for( int j = 0; j < 5; j++ )
        {
            double ep = -1.0 * ( inter[j ] - value )*( inter[j ] - value ) / sigma;
            weights[ j ] = exp( ep );
            
            sum += weights[ j ]*weights[ j ];
        }
        
        for( int j = 0; j < 5; j++ )
        {
            weights[ j ] =( weights[j] * weights[j] ) / sum;
        }
        
        colors[ i ] = weights[ 0 ] * red + weights[ 1 ] * yellow + weights[ 2 ]*green + weights[ 3 ] *cyan + weights[ 4 ]*blue;

    }
}





