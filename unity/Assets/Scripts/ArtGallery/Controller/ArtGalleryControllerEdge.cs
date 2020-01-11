namespace ArtGallery
{

    using System.Collections.Generic;
    using UnityEngine;
    using Util.Geometry;
    using Util.Geometry.DCEL;
    using Util.Geometry.Polygon;
    using Main;
    using System.Linq;
    using System;
    using Util.Math;
    using Util.Algorithms.Polygon;

    public class ArtGalleryControllerEdge : ArtGalleryController
    {

        public Dictionary<LineSegment, ArtGalleryLightHouse> segmentsWithLighthouse = new Dictionary<LineSegment, ArtGalleryLightHouse>();
        private Dictionary<int, Face> faceIDs = new Dictionary<int, Face>();
        private Dictionary<LineSegment, int> edgeIDs = new Dictionary<LineSegment, int>();
        private Dictionary<int, HashSet<int>> visibleCompIDsPerEdgeID = new Dictionary<int, HashSet<int>>();

        //Unity references
        private VisibilityAreaDrawer m_areaDrawer = null;

        // Update is called once per frame COULD GIVE PROBLEMS - TEST NECESSARY
        void Update()
        {
            //handle input key presses
            if (m_areaDrawer != null)
            {
                if (Input.GetKeyDown("a"))
                {
                    m_areaDrawer.ToggleDrawAll();
                }
                if (Input.GetKeyDown("e"))
                {
                    m_areaDrawer.ToggleDrawEdges();
                }
                if (Input.GetKeyDown("v"))
                {
                    m_areaDrawer.ToggleDrawVertices();
                }
                if (Input.GetKeyDown("f"))
                {
                    m_areaDrawer.ToggleDrawFaces();
                }
            }
        }

        public override void InitLevel()
        {
            base.InitLevel();
            Debug.Log("Initialising Level Drawer...");
            m_areaDrawer = FindObjectOfType<VisibilityAreaDrawer>();

            //DCEL dcell = new DCEL();
            //foreach (LineSegment s in LevelPolygon.Segments) {
            //    dcell.AddSegment(s);
            //}
            DCEL dcell = computeVisibilityRegions(LevelPolygon);
            //Debug.Log(dcell);
            if (m_areaDrawer != null)
            {
                m_areaDrawer.VisibilityAreas = dcell;
            }
            Debug.Log("Level Drawer Initalised!");

            int five = calcNeededNrOfEdgeGuards();
        }

        /// <summary>
        /// Computes the visibilty regions of a given polygon and outputs it as a DCEL
        /// </summary>
        /// <param name="poly">The polygon to compute the visibility regions for</param>
        /// <returns>
        /// A DCEL containing the visibility regions of a polygon
        /// </returns>
        private DCEL computeVisibilityRegions(Polygon2D poly)
        {
            ICollection<LineSegment> segments = getVisibilitySegments(poly);

            Debug.Log("Number of vertices: " + LevelPolygon.VertexCount);
            Debug.Log("Number of vertices in Vertices: " + LevelPolygon.Vertices.Count);

            //foreach (Vector2 v in LevelPolygon.Vertices) {
            //    Debug.Log(v);
            //}

            DCEL temp = new DCEL();
            foreach (LineSegment s in segments)
            {
                temp.AddSegment(s);
                //Debug.Log(s);
            }

            return temp; //mergeSegments(segments);
        }

        /// <summary>
        /// Gets the (overlapping) Line Segments that are needed to compute a given polygon's visibility regions
        /// </summary>
        /// <param name="poly">The polygon to compute the segments for</param>
        /// <returns>
        /// The (overlapping) Line Segments that are needed to compute the given polygon's visbility regions
        /// </returns>
        private ICollection<LineSegment> getVisibilitySegments(Polygon2D poly)
        {
            LinkedList<LineSegment> segments = new LinkedList<LineSegment>();

            foreach (Vector2 v in poly.Vertices)
            {
                foreach (Vector2 v2 in poly.Vertices)
                {
                    if (v.Equals(v2))
                    {
                        continue;
                    }

                    if (poly.isConvex(v2) == false)
                    { // v2 is Reflex

                        LineSegment s = new LineSegment(v, v2);
                        if (isIntersectPolygon(poly, s))
                        {
                            // Skip if the two vertices cannot see each other.
                            // NOTE: This includes colinear points (meaning that only 'adjacent' colinear points are not skipped)
                            continue;
                        }

                        Vector2? closestIntersection = intersectPolygonClosest(poly, new Line(v, v2));

                        if (closestIntersection != null)
                        {
                            segments.AddFirst(new LineSegment(v2, (Vector2)closestIntersection));
                        }

                    }
                }
            }
            return segments;
        }

        /// <summary>
        /// Finds the closest intersection of a Line's 'second point' and a polygon.
        /// </summary>
        /// <param name="poly">The polygon to find the closest intersection with</param>
        /// <param name="l">A line segment to check intersections with. Only the Ray from the line's 'second' point is 
        /// considered for intersections. 
        /// The 'second point' of l is a concave vertex.
        /// TODO: change to Ray2D rather than Line</param>
        /// <returns>
        /// The closest intersection from the line's second point with the given polygon
        /// </returns>
        // Finds the closest proper intersection
        // l.p1 is vertex
        // l.p2 is concave vertex
        private Vector2? intersectPolygonClosest(Polygon2D poly, Line l)
        {
            LineSegment smallestSegment = null;
            foreach (LineSegment s in poly.Segments)
            {
                // Will contain the intersection of l and s (if any)
                Vector2? intersection = null;

                // Ignore segments that the ray starts in (if any)
                if (s.IsOnSegment(l.Point2))
                {
                    // Exception: if a segment starts in l2 and is on the line, there is no intersection
                    if (Line.Colinear(l.Point1, s.Point1, s.Point2)
                            && !(l.Point1 == s.Point1 || l.Point1 == s.Point2))
                    {
                        return null;
                    }

                    continue;
                }

                // Edge Case: Handle segments that completely overlap the given line (I.e. there are 
                // infinitely many intersections)
                if (l.IsOnLine(s.Point1) && l.IsOnLine(s.Point2))
                {
                    Debug.Log("xyzSegment overlaps line (" + l + "): " + s.Point1 + " | " + s.Point2);
                    // All points of the segment overlap!
                    if ((new LineSegment(l.Point2, s.Point1)).IsOnSegment(s.Point2))
                    {
                        // Endpoint of the segment is closer
                        intersection = s.Point2;
                    }
                    else if ((new LineSegment(l.Point2, s.Point2)).IsOnSegment(s.Point1))
                    {
                        // Beginpoint of the segment is closer
                        intersection = s.Point1;
                    }
                    //} else if (s.IsOnSegment(l.Point2)) {
                    //    continue;
                }
                else
                {
                    // A single point of the segment can overlap
                    intersection = s.Intersect(l);
                }

                // Ignore intersections on the 'wrong' side
                // TODO: optimize this(?); ignore segments on one side of the line perpendicular to l
                if (intersection != null && (new LineSegment(l.Point2, (Vector2)intersection)).IsOnSegment(l.Point1))
                {
                    continue;
                }

                // There is an intersection, and it was closer than the previous closest intersection
                if (intersection != null && (smallestSegment == null || smallestSegment.IsOnSegment((Vector2)intersection)))
                {
                    smallestSegment = new LineSegment(l.Point2, (Vector2)intersection);
                }
            }

            if (smallestSegment != null && !smallestSegment.Point2.Equals(l.Point2))
            {
                return smallestSegment.Point2;
            }
            return null;
        }


        /// <summary>
        /// Checks whether a given segment between a polygon's vertices intersects at least one of that polygon's segments
        /// </summary>
        /// <param name="poly">The polygon to check for intersections with</param>
        /// <param name="vertexSegment">A line segment from one of poly's vertices to another of poly's vertices
        /// (possibly overlapping one or more segments of poly)</param>
        /// <returns>
        /// true iff the given segment overlaps at least one of the polygon's segments of vertices,
        /// ignoring any segments with a shared begin or end point with vertexSegment
        /// </returns>
        private bool isIntersectPolygon(Polygon2D poly, LineSegment vertexSegment)
        {
            foreach (LineSegment s in poly.Segments)
            {
                // Skip the segments that have a begin or end point in common with the given segment
                if (s.Point1 == vertexSegment.Point1 || s.Point2 == vertexSegment.Point2
                        || s.Point1 == vertexSegment.Point2 || s.Point2 == vertexSegment.Point1)
                {
                    continue;
                }

                // Edge Case: Skip segments whose begin and end points are BOTH colinear with the given segment
                //if (Line.Colinear(s.Point1, s.Point2, vertexSegment.Point1)
                //        && Line.Colinear(s.Point1, s.Point2, vertexSegment.Point2)) {
                //    continue;
                //}

                // Skip begin and end points of segments
                Vector2? x = vertexSegment.Intersect(s);
                if (x != null)
                {
                    return true;
                }
            }
            return false;
        }

        /// <summary>
        /// Creates a DCEL from a given list of Line Segments, where intersections of those
        /// segments are transformed into DCEL vertices and where line segments between those
        /// intersections are transformed to DCEL edges.
        /// </summary>
        /// <param name="segments">A list of line segments</param>
        /// <returns>
        /// A DCEL created from the given list of segments
        /// </returns>
        private DCEL mergeSegments(ICollection<LineSegment> segments)
        {

            List<LineSegment> lineSegments = new List<LineSegment>();
            // For each pair of segments
            foreach (LineSegment s1 in segments)
            {
                foreach (LineSegment s2 in segments)
                {
                    // Segments must not be the same
                    if (!s1.Equals(s2))
                    {
                        Vector2? intersection = s1.Intersect(s2);

                        // TODO: Handle multiple intersections in the same linesegment

                        // If the segments intersect, create a new vertex at the intersecion and split the segments
                        if (intersection != null)
                        {
                            lineSegments.Add(new LineSegment(s1.Point1, (Vector2)intersection));
                            lineSegments.Add(new LineSegment(s1.Point2, (Vector2)intersection));
                            lineSegments.Add(new LineSegment(s2.Point1, (Vector2)intersection));
                            lineSegments.Add(new LineSegment(s2.Point2, (Vector2)intersection));
                        }
                    }
                }
            }

            // Create the DCEL
            DCEL dcel = new DCEL();

            // Add the segments to the DCEL
            foreach (LineSegment segment in lineSegments)
            {
                dcel.AddSegment(segment);
            }
            return dcel;
        }

        /// <summary>
        /// Get a face by its ID
        /// </summary>
        /// <param name="id">An ID</param>
        /// <returns>
        /// The face corresponding to the given id
        /// </returns>
        private Face getFaceByID(int id)
        {
            Face face = faceIDs[id];
            return face;
        }

        /// <summary>
        /// Set IDs for a list of faces
        /// </summary>
        /// <param name="faces">A list of faces</param>
        private void setFaceIDs(List<Face> faces)
        {
            for (int i = faceIDs.Count; i < faceIDs.Count + faces.Count - 1; i++)
            {
                faceIDs.Add(i, faces[0]);
                faces.RemoveAt(0);
            }
        }

        /// <summary>
        /// Get an edge by its ID
        /// </summary>
        /// <param name="id">An ID</param>
        /// <returns>
        /// The edge corresponding to the given id
        /// </returns>
        private int getEdgeIDByEdge(LineSegment segment)
        {
            int edgeID = edgeIDs[segment];
            return edgeID;
        }

        /// <summary>
        /// Set IDs for a list of edges
        /// </summary>
        /// <param name="edges">A list of edges</param>
        private void setEdgeIDs(List<LineSegment> edges)
        {
            Debug.Log("size of the list: " + edges.Count);
            int count = 0;
            int edgeIDCount = edgeIDs.Count;
            for (int i = edgeIDCount; i < edgeIDCount + edges.Count; i++)
            {
                Debug.Log(count);
                LineSegment currentEdge = edges[count];
                edgeIDs.Add(currentEdge, i);
                count++;
            }
        }



        // Source: http://tripsintech.com/rotate-a-point-around-another-point-code/
        /// <summary>
        /// Rotates 'p1' about 'p2' by 'angle' degrees clockwise.
        /// </summary>
        /// <param name="p1">Point to be rotated</param>
        /// <param name="p2">Point to rotate around</param>
        /// <param name="angle">Angle in degrees to rotate clockwise</param>
        /// <returns>The rotated point</returns>
        public Vector2 RotatePoint(Vector2 p1, Vector2 p2, double angle)
        {

            double radians = ConvertToRadians(angle);
            double sin = Math.Sin(radians);
            double cos = Math.Cos(radians);

            // Translate point back to origin
            p1.x -= p2.x;
            p1.y -= p2.y;

            // Rotate point
            double xnew = p1.x * cos - p1.y * sin;
            double ynew = p1.x * sin + p1.y * cos;

            // Translate point back
            Vector2 newPoint = new Vector2((int)xnew + p2.x, (int)ynew + p2.y);
            return newPoint;
        }

        // Source: http://tripsintech.com/rotate-a-point-around-another-point-code/
        /// <summary>
        /// Convert angle to radians
        /// </summary>
        /// <param name="angle">The angle to convert to radians</param>
        /// <returns>Returns the angle converted to radians</returns>
        public double ConvertToRadians(double angle)
        {
            return (Math.PI / 180) * angle;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="visiblePoints"></param>
        /// <param name="vertices"></param>
        /// <param name="nextElement"></param>
        /// <param name="z"></param>
        /// <param name="xAxis"></param>
        /// <returns></returns>
        private int newEdgeCrossesXAxis(Stack<Vector2> visiblePoints, List<Vector2> vertices, int nextElement, Vector2 z, Line xAxis)
        {
            // Whenever the stack content changes from (s_0, s_1, ..., s_m) to (s_0, s_1, ..., s_(m+1)) such that
            // line segment (s_m, s_(m+1)) crosses the x-axis at v
            LineSegment segmentS = new LineSegment(visiblePoints.ElementAt(1), visiblePoints.Peek());
            Vector2? v = segmentS.Intersect(xAxis);
            if (v != null)
            {
                // Continue scanning vertices of P until an edge (u_k, u_(k+1)) intersecting line segment (s_0, v) is found
                // Then the treatment for this situation is exactly the same as region3. 
                for (int i = 0; i < vertices.Count; i++)
                {
                    int index = (nextElement + i) % vertices.Count;
                    int indexNext = (nextElement + i - 1) % vertices.Count;
                    LineSegment segmentK = new LineSegment(vertices[index], vertices[indexNext]);
                    if (segmentK.Intersect(new LineSegment(visiblePoints.Last(), (Vector2)v)) != null) ;
                    {
                        return index;
                    }
                }
            }
            return -1;
        }

        /// <summary>
        /// Handles region 1 in case 1 of algorithm VP(P, z)
        /// </summary>
        /// <param name="visiblePoints">Current stack of "visible" points</param>
        /// <param name="vertices">List of vertices in CCW order</param>
        /// <param name="nextElement">Next element to check</param>
        /// <param name="z">Point of origin</param>
        /// <param name="xAxis">Horizontal line through z</param>
        /// <returns>Returns a stack of visible points in P from point z</returns>
        private Stack<Vector2> case1Region1(Stack<Vector2> visiblePoints, List<Vector2> vertices, int nextElement, Vector2 z, Line xAxis)
        {
            // NextElement in Region1
            // NextElement not visible from z
            // Scan vertices of P from nextElement until some vertex u_k
            int crossedXAxisCount = 0;
            for (int i = 0; i < vertices.Count; i++)
            {
                int index = (nextElement + i) % (vertices.Count);
                int indexPrev = (nextElement + i - 1) % (vertices.Count);
                if (indexPrev < 0) // indexPrev can only be negative for -1. Then the previous index is the index of the last vertex in the list
                {
                    indexPrev = vertices.Count - 1;
                }
                LineSegment segment = new LineSegment(vertices[indexPrev], vertices[index]);
                
                Ray2D halfLine = new Ray2D(visiblePoints.Peek(), RotatePoint(z, visiblePoints.Peek(), 180));
                Vector2? v = segment.Intersect(halfLine);
                //If the x-axis intersection count is even and the intersection is not the origin of the ray
                if (!v.Equals(null) && !MathUtil.EqualsEps((Vector2)v, visiblePoints.Peek()) && crossedXAxisCount % 2 == 0)
                {
                    nextElement = (index + 1) % vertices.Count;
                    visiblePoints.Push((Vector2)v);
                    if (index == 0) // Algorithm terminates when u_0 is pushed on the stack again 
                    {
                        return visiblePoints;
                    }
                    visiblePoints.Push(vertices[index]);
                    int newIndex = newEdgeCrossesXAxis(visiblePoints, vertices, nextElement, z, xAxis);
                    if (newIndex >= 0)
                    {
                        nextElement = newIndex;
                        return region3(visiblePoints, vertices, nextElement, z, xAxis);
                    }
                    break;
                }
                // TODO check whether positive x-axis is correct
                // If segment intersects the positive x-axis then we increase the count
                if (segment.Intersect(xAxis)?.x > z.x)
                {
                    crossedXAxisCount++;
                }
            }
            return case1(visiblePoints, vertices, nextElement, z, xAxis);
        }

        /// <summary>
        /// Handles region 2 in case 1 of algorithm VP(P, z)
        /// </summary>
        /// <param name="visiblePoints">Current stack of "visible" points</param>
        /// <param name="vertices">List of vertices in CCW order</param>
        /// <param name="nextElement">Next element to check</param>
        /// <param name="z">Point of origin</param>
        /// <param name="xAxis">Horizontal line through z</param>
        /// <returns>Returns a stack of visible points in P from point z</returns>
        private Stack<Vector2> case1Region2(Stack<Vector2> visiblePoints, List<Vector2> vertices, int nextElement, Vector2 z, Line xAxis)
        {
            // NextElement in Region2
            if (nextElement == 0) // Algorithm terminates when u_0 is pushed on the stack again
            {
                return visiblePoints;
            }
            visiblePoints.Push(vertices[nextElement]);
            nextElement = (nextElement + 1) % vertices.Count;
            int newIndex = newEdgeCrossesXAxis(visiblePoints, vertices, nextElement, z, xAxis);
            if (newIndex >= 0)
            {
                nextElement = newIndex;
                return region3(visiblePoints, vertices, nextElement, z, xAxis);
            }
            return case1(visiblePoints, vertices, nextElement, z, xAxis);
        }

        /// <summary>
        /// Handles region 1 in case 2 of algorithm VP(P, z)
        /// </summary>
        /// <param name="visiblePoints">Current stack of "visible" points</param>
        /// <param name="vertices">List of vertices in CCW order</param>
        /// <param name="nextElement">Next element to check</param>
        /// <param name="z">Point of origin</param>
        /// <param name="xAxis">Horizontal line through z</param>
        /// <returns>Returns a stack of visible points in P from point z</returns>
        private Stack<Vector2> case2Region1(Stack<Vector2> visiblePoints, List<Vector2> vertices, int nextElement, Vector2 z, Line xAxis)
        {
            // Scanning vertices of P until vertex u_k such that line segment (u_(k-1), u_k) intersects line segment (s_(j-1), s_j) at v for the first time
            for (int i = 0; i < vertices.Count; i++)
            {
                int index = (nextElement + i) % vertices.Count;
                int indexPrev = (nextElement + i - 1) % vertices.Count;
                if (indexPrev < 0) // indexPrev can only be negative for -1. Then the previous index is the index of the last vertex in the list
                {
                    indexPrev = vertices.Count - 1;
                }
                LineSegment segmentK = new LineSegment(vertices[indexPrev], vertices[index]);
                Vector2? v = segmentK.Intersect(new LineSegment(visiblePoints.ElementAt(1), visiblePoints.Peek()));
                if (v != null)
                {
                    // Remaining steps identical to region 3 where s_(j-1) is v
                    Vector2 topOfStack = visiblePoints.Pop(); 
                    visiblePoints.Pop(); // We remove the element previously at s_(j-1). otherwise C5 in the example would contain v2
                    visiblePoints.Push((Vector2)v); 
                    visiblePoints.Push(topOfStack);
                    nextElement = index; // We update the index, otherwise v_3 would be in C5 instead of v_4
                    break;
                }
            }
            return region3(visiblePoints, vertices, nextElement, z, xAxis);
        }

        /// <summary>
        /// Handles region 2 in case 2 of algorithm VP(P, z)
        /// </summary>
        /// <param name="visiblePoints">Current stack of "visible" points</param>
        /// <param name="vertices">List of vertices in CCW order</param>
        /// <param name="nextElement">Next element to check</param>
        /// <param name="z">Point of origin</param>
        /// <param name="xAxis">Horizontal line through z</param>
        /// <returns>Returns a stack of visible points in P from point z</returns>
        private Stack<Vector2> case2Region2(Stack<Vector2> visiblePoints, List<Vector2> vertices, int nextElement, Vector2 z, Line xAxis)
        {
            // The configuration becomes (u_(i+1); s_0, s_1, ..., s_j, u_i) and belongs to C_1
            if (nextElement == 0) // Algorithm terminates when u_0 is pushed on the stack again
            {
                return visiblePoints;
            }
            visiblePoints.Push(vertices[nextElement]);
            nextElement = (nextElement + 1) % vertices.Count;
            int newIndex = newEdgeCrossesXAxis(visiblePoints, vertices, nextElement, z, xAxis);
            if (newIndex >= 0)
            {
                nextElement = newIndex;
                return region3(visiblePoints, vertices, nextElement, z, xAxis);
            }
            return case1(visiblePoints, vertices, nextElement, z, xAxis);
        }

        /// <summary>
        /// Handles region 3 of cases 1 and 2 of algorithm VP(P, z)
        /// </summary>
        /// <param name="visiblePoints">Current stack of "visible" points</param>
        /// <param name="vertices">List of vertices in CCW order</param>
        /// <param name="nextElement">Next element to check</param>
        /// <param name="z">Point of origin</param>
        /// <param name="xAxis">Horizontal line through z</param>
        /// <returns>Returns a stack of visible points in P from point z</returns>
        private Stack<Vector2> region3(Stack<Vector2> visiblePoints, List<Vector2> vertices, int nextElement, Vector2 z, Line xAxis)
        {
            // NextElement in Region3
            // Edge (NextElement-1, NextElement) blocks points in the stack
            // Pop elements off the stack until some s_m such that edge (u_(i-1), u_i)
            // intersects (s_m, s_(m+1)) at v or u_i lies to the left of line (z, s_m)
            LineSegment currNextLine = new LineSegment(vertices[nextElement - 1], vertices[nextElement]);
            Vector2? v = currNextLine.Intersect(new LineSegment(visiblePoints.Peek(), visiblePoints.ElementAt(1)));
            while (v == null || !new Line(z, visiblePoints.ElementAt(1)).PointRightOfLine(vertices[nextElement]))
            {
                visiblePoints.Pop();
                v = currNextLine.Intersect(new LineSegment(visiblePoints.Peek(), visiblePoints.ElementAt(1)));
            }
            // Still have s_(m+1) on the stack, which needs to be removed from the stack
            Vector2 sm1 = visiblePoints.Pop();
            // If edge (u_(i-1), u_i) intersects (s_m, s_(m+1))
            if (v != null)
            {
                LineSegment edgeSMV = new LineSegment(visiblePoints.Peek(), (Vector2)v); //(s_m, v)
                // Scan vertices of P in CCW order until edge (u_(k-1), u_k) crosses edge (s_m, v) at w
                // for the first time. 
                for (int i = 0; i < vertices.Count; i++)
                {
                    int index = (nextElement + i) % (vertices.Count); //k
                    int indexPrev = (nextElement + i - 1) % (vertices.Count); //k-1
                    if (indexPrev < 0) // indexPrev can only be negative for -1. Then the previous index is the index of the last vertex in the list
                    {
                        indexPrev = vertices.Count - 1;
                    }
                    LineSegment segment = new LineSegment(vertices[indexPrev], vertices[index]); //(u_(k-1), u_k)
                    Vector2? w = segment.Intersect(edgeSMV);
                    if (w != null)
                    {
                        // New configuration becomes:
                        // (u_(k+1) ; s_0, s_1, ..., s_m, w, u_k)
                        nextElement = (index + 1) % (vertices.Count);
                        visiblePoints.Push((Vector2)w);
                        if (index == 0) // Algorithm terminates when u_0 is pushed on the stack again
                        {
                            return visiblePoints;
                        }
                        visiblePoints.Push(vertices[index]);
                        int newIndex = newEdgeCrossesXAxis(visiblePoints, vertices, nextElement, z, xAxis);
                        if (newIndex >= 0)
                        {
                            nextElement = newIndex;
                            return region3(visiblePoints, vertices, nextElement, z, xAxis);
                        }
                        break;
                    }
                }
                return case1(visiblePoints, vertices, nextElement, z, xAxis);
            }
            else // u_i lies to the left of line (z, s_m)
            {
                v = new LineSegment(visiblePoints.Peek(), sm1).Intersect(new Line(z, vertices[nextElement]));
                visiblePoints.Push((Vector2)v);
                if (nextElement == 0) // Algorithm terminates when u_0 is pushed on the stack again
                {
                    return visiblePoints;
                }
                visiblePoints.Push(vertices[nextElement]);
                int newIndex = newEdgeCrossesXAxis(visiblePoints, vertices, nextElement, z, xAxis);
                if (newIndex >= 0)
                {
                    nextElement = newIndex;
                    return region3(visiblePoints, vertices, nextElement, z, xAxis);
                }
                nextElement = (nextElement + 1) % vertices.Count;
                return case2(visiblePoints, vertices, nextElement, z, xAxis);
            }
        }

        /// <summary>
        /// Handles case 1 of algorithm VP(P, z)
        /// </summary>
        /// <param name="visiblePoints">Current stack of "visible" points</param>
        /// <param name="vertices">List of vertices in CCW order</param>
        /// <param name="nextElement">Next element to check</param>
        /// <param name="z">Point of origin</param>
        /// <param name="xAxis">Horizontal line through z</param>
        /// <returns>Returns a stack of visible points in P from point z</returns>
        private Stack<Vector2> case1(Stack<Vector2> visiblePoints, List<Vector2> vertices, int nextElement, Vector2 z, Line xAxis)
        {
            // Polygon from z and the current stack
            Polygon2D polyZAndStack = new Polygon2D(visiblePoints);
            polyZAndStack.AddVertexFirst(z);
            // Determine region in which nextElement lies
            // Region3 is the interior of Z
            // Region1 and Region2 are created from the exterior of Z
            // by a line through z and the first element on the stack
            if (polyZAndStack.ContainsInside(vertices[nextElement]))
            {
                return region3(visiblePoints, vertices, nextElement, z, xAxis);
            }
            else
            {
                // NextElement in Region1 or Region2
                Line regionSplit = new Line(z, visiblePoints.Peek());
                if (regionSplit.PointRightOfLine(vertices[nextElement]))
                {
                    return case1Region1(visiblePoints, vertices, nextElement, z, xAxis);
                }
                else
                {
                    return case1Region2(visiblePoints, vertices, nextElement, z, xAxis);
                }
            }
        }

        /// <summary>
        /// Computes the "chain" of vertices starting from a startvertex up to and including an endvertex
        /// </summary>
        /// <param name="vertices">List of vertices from a polygon in counterclockwise order</param>
        /// <param name="startVertex">The startvertex</param>
        /// <param name="endVertex">The endvertex</param>
        /// <returns>Returns the "chain" of vertices starting from a startvertex up to and including an endvertex</returns>
        private List<Vector2> chain(List<Vector2> vertices, Vector2 startVertex, Vector2 endVertex)
        {
            List<Vector2> result = new List<Vector2>();
            int startIndex = vertices.IndexOf(startVertex);
            int endIndex = vertices.IndexOf(endVertex);
            for (int i = 0; i < vertices.Count; i++)
            {
                int index = (startIndex + i) % vertices.Count;
                if (index == endIndex)
                {
                    result.Add(vertices[index]);
                    break;
                }
                else
                {
                    result.Add(vertices[index]);
                }
            }
            return result;
        }

        /// <summary>
        /// Handles case 2 of algorithm VP(P, z)
        /// </summary>
        /// <param name="visiblePoints">Current stack of "visible" points</param>
        /// <param name="vertices">List of vertices in CCW order</param>
        /// <param name="nextElement">Next element to check</param>
        /// <param name="z">Point of origin</param>
        /// <param name="xAxis">Horizontal line through z</param>
        /// <returns>Returns a stack of visible points in P from point z</returns>
        private Stack<Vector2> case2(Stack<Vector2> visiblePoints, List<Vector2> vertices, int nextElement, Vector2 z, Line xAxis)
        {
            // Three regions: 
            // R1, defined by Ch(s_(j-1), s_j) (chain from s_(j-1) to s_j) and line segment (s_j, s_(j-1))
            // R3, which is the interior of the polygon defined by line segment (z, s_0), followed by Ch(s_0, s_(j-1)), and line segment (s_(j-1), z)
            // R2, which is the remaining region

            // Defining R1
            // TODO FIX THESE POLYGONS, CURRENTLY WRONG
            List<Vector2> chainR1 = chain(vertices, visiblePoints.ElementAt(1), visiblePoints.Peek());
            Polygon2D polyR1 = new Polygon2D(chainR1);
            polyR1.AddVertexAfter(visiblePoints.Peek(), chainR1.Last());
            polyR1.AddVertexAfter(visiblePoints.ElementAt(1), visiblePoints.Peek());

            // Defining R3
            List<Vector2> chainR3 = chain(vertices, visiblePoints.Last(), visiblePoints.ElementAt(1));
            Polygon2D polyR3 = new Polygon2D();
            polyR3.AddVertexFirst(z);
            polyR3.AddVertex(visiblePoints.Last());
            foreach (var vertex in chainR3)
            {
                polyR3.AddVertex(vertex);
            }

            if (polyR1.ContainsInside(vertices[nextElement]))
            {
                // R1
                return case2Region1(visiblePoints, vertices, nextElement, z, xAxis);
            }
            else if (polyR3.ContainsInside(vertices[nextElement]))
            {
                // R3
                // Treatment for this case identical to treatment for region 3 of case 1
                return region3(visiblePoints, vertices, nextElement, z, xAxis);
            }
            else
            {
                // R2
                return case2Region2(visiblePoints, vertices, nextElement, z, xAxis);
            }
        }

        /// <summary>
        /// Computes the weakly visible points in polygon P from point z
        /// Source: https://www.sciencedirect.com/science/article/pii/0734189X83900658
        /// </summary>
        /// <param name="polygon">A polygon</param>
        /// <param name="z">A point</param>
        /// <returns>A list of visible points in polygon P from point z</returns>
        private List<Vector2> computeWeaklyVisiblePointsInPFromZ(Polygon2D polygon, Vector2 z)
        {
            // Stack containing the vertices and some boundary points of P that are "visible" from z so far and at 
            // termination contains all visible vertices from z
            Stack<Vector2> visiblePoints = new Stack<Vector2>();
            // Let z be the origin
            Line xAxis = new Line(z, 0f);
            // Let u0 be the intersection of an edge of P and the positive x axis that has the smallest x coordinate
            Vector2? u0 = null;
            Vector2? nextVertex = null;
            foreach (var segment in polygon.Segments)
            {
                var intersection = segment.Intersect(xAxis);
                if (intersection != null)
                {
                    if (u0 == null && intersection?.x > z.x)
                    {
                        u0 = intersection;
                        nextVertex = segment.Point2;
                    }
                    else
                    {
                        if (intersection?.x > z.x && intersection?.x < u0?.x)
                        {
                            u0 = intersection;
                            nextVertex = segment.Point2;
                        }
                    }
                }
            }

            // Let u0 be a vertex of polygon P, and the vertices are ordered counterclockwise
            List<Vector2> vertices = new List<Vector2>(polygon.Vertices.ToList());
            int nextIndex = vertices.IndexOf((Vector2)nextVertex);
            vertices.Insert(nextIndex, (Vector2)u0);
            // Shift the list such that u0 is at the end of the list
            while (nextIndex >= 0)
            {
                var first = vertices[0];
                vertices.Remove(first);
                vertices.Add(first);
                nextIndex--;
            }
            // Reverse the list such that the vertices are ordered counter clockwise
            // Note that u0 is now at the front of the list
            vertices.Reverse();

            //DEBUG CODE 
            try
            {
                if (MathUtil.EqualsEps((Vector2)u0, vertices[0]))
                {
                    throw new Exception("Exception: u0 is not at the front of the list of vertices");
                }
            }
            catch (Exception e)
            {
                Debug.LogException(e, this);
            }
            // DEBUG CODE 


            // Initially the stack contains u0 and u1
            visiblePoints.Push((Vector2)u0);
            visiblePoints.Push(vertices[1]);
            // Next element to check is u2
            int nextElement = vertices.IndexOf(vertices[2]);

            // Get index of element u_{i-2}
            int element2Prev;
            if (nextElement - 2 < 0)
            {
                element2Prev = vertices.Count - (nextElement - 2);
            }
            else
            {
                element2Prev = nextElement - 2;
            }

            //Check if u_{i-2} lies to the right or left of line segment (z, s_j)
            Stack<Vector2> result;
            if (new LineSegment(z, visiblePoints.Peek()).IsRightOf(vertices[element2Prev]))
            {
                // Case C1
                result = case1(visiblePoints, vertices, nextElement, z, xAxis);
            }
            else
            {
                // Case C2
                result = case2(visiblePoints, vertices, nextElement, z, xAxis);
            }
            return result.ToList();
        }

        /// <summary>
        /// Compute the visible components per edge
        /// </summary>
        private void computeVisibleComponentsPerEdge(DCEL dcel)
        {
            foreach (var edgeID in edgeIDs)
            {
                // Store visible components(/faces) for this edge by index
                HashSet<int> visibleComps = new HashSet<int>();
                // Compute visible components (faces)
                foreach (var faceID in faceIDs)
                {
                    // For every convex component c, take one of its vertices z
                    DCELVertex vertex = faceID.Value.OuterVertices.ToList()[0];
                    // Compute VP(P, z); weakly visible points of P from z
                    List<Vector2> weakVisiblePoints = computeWeaklyVisiblePointsInPFromZ(LevelPolygon, vertex.Pos);
                    // For every vertex v in VP(P, z), add c to the set of the edge containing v
                    foreach (var v in weakVisiblePoints)
                    {
                        foreach (var segment in dcel.Edges)
                        {
                            if (segment.Segment.IsOnSegment(v))
                            {
                                visibleComps.Add(faceID.Key);
                            }
                        }
                    }
                }
                visibleCompIDsPerEdgeID.Add(edgeID.Value, visibleComps);
            }
        }

        /// <summary>
        /// Approximate the minimum number of edge guards needed for the current level's polygon
        /// with approximation ratio O(log(n))
        /// </summary>
        /// <returns>
        /// An approximation of the minimum number of edge guards needed for the current level's 
        /// polygon with approximation ratio O(log(n))
        /// </returns>
        private int calcNeededNrOfEdgeGuards()
        {
            //// Draw lines through every pair of vertices 
            //List<Line> lines = generateLines(LevelPolygon.Vertices.ToList());
            //// Create convex components 
            //List<LineSegment> lineSegments = generateLineSegments(lines);
            //DCEL dcel = createDCELFromPolygonAndSegments(LevelPolygon, lineSegments);
            //List<Face> faces = dcel.InnerFaces.ToList();
            //// Give IDs to the faces
            //setFaceIDs(faces);
            //// Give IDs to the edges
            setEdgeIDs(LevelPolygon.Segments.ToList());
            //// Store visible convex components(/faces) for each edge by index
            //computeVisibleComponentsPerEdge(dcel);
            //// Calculate the needed number of edge guards
            //return SetCover.Solve(new HashSet<int>(faceIDs.Keys), visibleCompIDsPerEdgeID);
            return 5;
        }

        public override void HandleIslandClick()
        {
            // obtain mouse position
            var worldlocation = Camera.main.ScreenPointToRay(Input.mousePosition).origin;
            worldlocation.z = -2f;

            //Calculate nearest segment from island click
            float minDistance = float.MaxValue;
            LineSegment closestSegment = null;

            var segments = LevelPolygon.Segments;
            foreach (var segment in segments)
            {
                Vector3 point1Vector3 = segment.Point1;
                point1Vector3.z = -2f;
                Vector3 point2Vector3 = segment.Point2;
                point2Vector3.z = -2f;
                float distanceToSegment = UnityEditor.HandleUtility.DistancePointLine(worldlocation, point1Vector3, point2Vector3);
                if (distanceToSegment < minDistance)
                {
                    minDistance = distanceToSegment;
                    closestSegment = segment;
                }
            }

            if (segmentsWithLighthouse.ContainsKey(closestSegment))
            {
                var lighthouseToRemove = segmentsWithLighthouse[closestSegment];

                // destroy the lighthouse
                m_solution.RemoveLighthouse(lighthouseToRemove);
                Destroy(lighthouseToRemove.gameObject);
                UpdateLighthouseText();

                segmentsWithLighthouse.Remove(closestSegment);
                CheckSolution();
                return;
            }

            // return if lighthouse was already selected or player can place no more lighthouses
            if (m_selectedLighthouse != null || m_solution.Count >= m_maxNumberOfLighthouses)
                return;

            Vector3 locationForLighthouse = closestSegment.Midpoint;
            locationForLighthouse.z = -2f;

            // create a new lighthouse from prefab
            var go = Instantiate(m_lighthousePrefab, locationForLighthouse, Quaternion.identity) as GameObject;

            // Add closest line segment to lighthouse
            go.GetComponent<ArtGalleryLightHouse>().m_segment = closestSegment;

            go.GetComponent<ArtGalleryLightHouse>().UpdateVision();

            // add lighthouse to art gallery solution
            m_solution.AddLighthouse(go);
            UpdateLighthouseText();

            // Add the segment/lighthouse combination to segmentWithLighthouse
            segmentsWithLighthouse.Add(closestSegment, go.GetComponent<ArtGalleryLightHouse>());

            CheckSolution();
        }

        /// <summary>
        /// Update the vision polygon for the given lighthouse.
        /// Calculates the visibility polygon.
        /// </summary>
        /// <param name="m_lighthouse"></param>
        public override void UpdateVision(ArtGalleryLightHouse m_lighthouse)
        {
            if (LevelPolygon.ContainsInside(m_lighthouse.Pos))
            {
                LineSegment selectedEdge = m_lighthouse.m_segment;
                int edgeID = getEdgeIDByEdge(selectedEdge);
                HashSet<int> visibleCompIDs = visibleCompIDsPerEdgeID[edgeID];
                List<Face> faces = new List<Face>();
                foreach (var id in visibleCompIDs)
                {
                    faces.Add(getFaceByID(id));
                }
                IEnumerable<Vector2> outerPointsAllFaces = new Vector2[] { };
                foreach (var face in faces)
                {
                    Enumerable.Concat(outerPointsAllFaces, face.OuterPoints);
                }
                Polygon2D vision = new Polygon2D(outerPointsAllFaces);
                // update lighthouse visibility
                m_lighthouse.VisionPoly = vision;
                m_lighthouse.VisionAreaMesh.Polygon = vision;
            }
            else
            {
                // remove visibility polygon from lighthouse
                m_lighthouse.VisionPoly = null;
                m_lighthouse.VisionAreaMesh.Polygon = null;
            }
        }
    }
}
