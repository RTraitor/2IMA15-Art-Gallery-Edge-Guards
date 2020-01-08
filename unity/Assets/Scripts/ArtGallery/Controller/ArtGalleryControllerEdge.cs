namespace ArtGallery {

    using System.Collections.Generic;
    using UnityEngine;
    using Util.Geometry;
    using Util.Geometry.DCEL;
    using Util.Geometry.Polygon;
    using Main;
    using System.Linq;
    using System;
    using Util.Math;

    public class ArtGalleryControllerEdge : ArtGalleryController
    {

        public Dictionary<LineSegment, ArtGalleryLightHouse> segmentsWithLighthouse = new Dictionary<LineSegment, ArtGalleryLightHouse>();
        private Dictionary<int, Face> faceIDs = new Dictionary<int, Face>();
        private Dictionary<int, LineSegment> edgeIDs = new Dictionary<int, LineSegment>();
        private Dictionary<int, HashSet<int>> visibleCompIDsPerEdgeID = new Dictionary<int, HashSet<int>>();

        // Update is called once per frame COULD GIVE PROBLEMS - TEST NECESSARY
        void Update() { }

        /// <summary>
        /// Generate lines through every pair of vertices from a pointset
        /// </summary>
        /// <param name="vertices">A list of vertices</param>
        /// <returns>
        /// A list of lines through each pair of vertices in the given pointset "vertices"
        /// </returns>
        private List<Line> generateLines(List<Vector2> vertices) 
        {
            List<Line> lines = new List<Line>();
            foreach (Vector2 vertex1 in vertices)
            {
                foreach (Vector2 vertex2 in vertices)
                {
                    if (vertex1 != vertex2)
                    {
                        Line line = new Line(vertex1, vertex2);
                        // Since each pair is checked twice, we check if we  
                        // already included the line instead of adding it twice
                        if (!lines.Contains(line))
                        {
                            lines.Add(line);
                        }
                    }
                }
            }
            return lines;
        }

        /// <summary>
        /// Line segments are created from the intersection of each pair of lines to each
        /// of the (at most four) points corresponding to the intersecting lines
        /// </summary>
        /// <param name="lines">A list of lines</param>
        /// <returns>
        /// Line segments created from the intersection of each pair of lines to each
        /// of the (at most four) points corresponding to the intersecting lines
        /// </returns>
        private List<LineSegment> generateLineSegments(List<Line> lines)
        {
            List<Vector2> intersections = new List<Vector2>();
            List<LineSegment> lineSegments = new List<LineSegment>();
            // For each pair of lines
            foreach (Line line1 in lines)
            {
                foreach (Line line2 in lines)
                {
                    if (!line1.Equals(line2))
                    {
                        Vector2 intersection = (Vector2)line1.Intersect(line2);
                        // Each pair of lines is checked twice, we check if we 
                        // already included the intersection instead of adding it twice
                        if (intersection != null && LevelPolygon.ContainsInside(intersection))
                        {
                            intersections.Add(intersection);
                            lineSegments.Add(new LineSegment(line1.Point1, intersection));
                            lineSegments.Add(new LineSegment(line1.Point2, intersection));
                            lineSegments.Add(new LineSegment(line2.Point1, intersection));
                            lineSegments.Add(new LineSegment(line2.Point2, intersection));
                        }
                    }
                }
            }
            return lineSegments;
        }

        /// <summary>
        /// Create a DCEL from a list of line segments and a polygon
        /// </summary>
        /// <param name="polygon">A Polygon</param>
        /// <param name="lineSegments">A list of line segments</param>
        /// <returns>
        /// A DCEL created from a starting polygon and a list of lineSegments
        /// </returns>
        private DCEL createDCELFromPolygonAndSegments(Polygon2D polygon, List<LineSegment> lineSegments)
        {
            DCEL dcel = new DCEL();
            // Add the segments of the polygon to the DCEL
            foreach (LineSegment segment in polygon.Segments)
            {
                dcel.AddSegment(segment);
            }

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
        private LineSegment getEdgeByID(int id)
        {
            LineSegment edge = edgeIDs[id];
            return edge;
        }

        /// <summary>
        /// Set IDs for a list of edges
        /// </summary>
        /// <param name="edges">A list of edges</param>
        private void setEdgeIDs(List<LineSegment> edges)
        {
            for (int i = edgeIDs.Count; i < edgeIDs.Count + edges.Count - 1; i++)
            {
                edgeIDs.Add(i, edges[i]);
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
        private Boolean newEdgeCrossesXAxis(Stack<Vector2> visiblePoints, List<Vector2> vertices, int nextElement, Vector2 z, Line xAxis)
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
                        return true;
                    }
                }
            }
            return false;
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
                LineSegment segment = new LineSegment(vertices[indexPrev], vertices[index]);
                // If segment intersects the positive x-axis then we increase the count
                if (segment.Intersect(xAxis)?.x > z.x)
                {
                    crossedXAxisCount++;
                }
                Ray2D halfLine = new Ray2D(visiblePoints.Peek(), RotatePoint(z, visiblePoints.Peek(), 180));
                Vector2? intersection = segment.Intersect(halfLine);
                //If the x-axis intersection count is even and the intersection is not the origin of the ray
                if (!intersection.Equals(null) && !MathUtil.EqualsEps((Vector2)intersection, visiblePoints.Peek()) && crossedXAxisCount % 2 == 0)
                {
                    Vector2 v = (Vector2)segment.Intersect(halfLine);
                    nextElement = (index + 1) % vertices.Count;
                    visiblePoints.Push(v);
                    if (newEdgeCrossesXAxis(visiblePoints, vertices, nextElement, z, xAxis))
                    {
                        return region3(visiblePoints, vertices, nextElement, z, xAxis);
                    }
                    if (index == 0) // Algorithm terminates when u_0 is pushed on the stack again
                    {
                        return visiblePoints;
                    }
                    visiblePoints.Push(vertices[index]);
                    if (newEdgeCrossesXAxis(visiblePoints, vertices, nextElement, z, xAxis))
                    {
                        return region3(visiblePoints, vertices, nextElement, z, xAxis);
                    }
                    break;
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
            if (newEdgeCrossesXAxis(visiblePoints, vertices, nextElement, z, xAxis))
            {
                return region3(visiblePoints, vertices, nextElement, z, xAxis);
            }
            nextElement = (nextElement + 1) % vertices.Count;
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
                LineSegment segmentK = new LineSegment(vertices[indexPrev], vertices[index]);
                Vector2? v = segmentK.Intersect(new LineSegment(visiblePoints.ElementAt(1), visiblePoints.Peek()));
                if (v != null)
                {
                    // Remaining steps identical to region 3 where s_(j-1) is v
                    Vector2 topOfStack = visiblePoints.Pop();
                    visiblePoints.Pop(); // TODO: This line might need to be removed, not sure
                    visiblePoints.Push((Vector2)v);
                    if (newEdgeCrossesXAxis(visiblePoints, vertices, nextElement, z, xAxis))
                    {
                        return region3(visiblePoints, vertices, nextElement, z, xAxis);
                    }
                    visiblePoints.Push(topOfStack);
                    if (newEdgeCrossesXAxis(visiblePoints, vertices, nextElement, z, xAxis))
                    {
                        return region3(visiblePoints, vertices, nextElement, z, xAxis);
                    }
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
            if (newEdgeCrossesXAxis(visiblePoints, vertices, nextElement, z, xAxis))
            {
                return region3(visiblePoints, vertices, nextElement, z, xAxis);
            }
            nextElement = (nextElement + 1) % vertices.Count;
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
                    LineSegment segment = new LineSegment(vertices[indexPrev], vertices[index]); //(u_(k-1), u_k)
                    Vector2? w = segment.Intersect(edgeSMV);
                    if (w != null)
                    {
                        // New configuration becomes:
                        // (u_(k+1) ; s_0, s_1, ..., s_m, w, u_k)
                        nextElement = (index + 1) % (vertices.Count);
                        visiblePoints.Push((Vector2)w);
                        if (newEdgeCrossesXAxis(visiblePoints, vertices, nextElement, z, xAxis))
                        {
                            return region3(visiblePoints, vertices, nextElement, z, xAxis);
                        }
                        if (index == 0) // Algorithm terminates when u_0 is pushed on the stack again
                        {
                            return visiblePoints;
                        }
                        visiblePoints.Push(vertices[index]);
                        if (newEdgeCrossesXAxis(visiblePoints, vertices, nextElement, z, xAxis))
                        {
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
                if (newEdgeCrossesXAxis(visiblePoints, vertices, nextElement, z, xAxis))
                {
                    return region3(visiblePoints, vertices, nextElement, z, xAxis);
                }
                if (nextElement == 0) // Algorithm terminates when u_0 is pushed on the stack again
                {
                    return visiblePoints;
                }
                visiblePoints.Push(vertices[nextElement]);
                if (newEdgeCrossesXAxis(visiblePoints, vertices, nextElement, z, xAxis))
                {
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
            List<Vector2> vertices = new List<Vector2>((List<Vector2>) polygon.Vertices);
            int nextIndex = vertices.IndexOf((Vector2) nextVertex);
            vertices.Insert(nextIndex, (Vector2) u0);
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

            // Initially the stack contains u0 and u1
            visiblePoints.Push((Vector2) u0);
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
                    DCELVertex vertex = ((List<DCELVertex>)faceID.Value.OuterVertices)[0];
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
                visibleCompIDsPerEdgeID.Add(edgeID.Key, visibleComps);
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
            // Draw lines through every pair of vertices 
            List<Line> lines = generateLines((List<Vector2>) LevelPolygon.Vertices);
            // Create convex components 
            List<LineSegment> lineSegments = generateLineSegments(lines);
            DCEL dcel = createDCELFromPolygonAndSegments(LevelPolygon, lineSegments);
            List<Face> faces = (List<Face>) dcel.InnerFaces;
            // Give IDs to the faces
            setFaceIDs(faces);
            // Give IDs to the edges
            setEdgeIDs((List<LineSegment>) LevelPolygon.Segments);
            // Store visible convex components(/faces) for each edge by index
            computeVisibleComponentsPerEdge(dcel);
            // Calculate the needed number of edge guards
            return SetCover.Solve(new HashSet<int>(faceIDs.Keys), visibleCompIDsPerEdgeID);
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

            // add lighthouse to art gallery solution
            m_solution.AddLighthouse(go);
            UpdateLighthouseText();

            // Add the segment/lighthouse combination to segmentWithLighthouse
            segmentsWithLighthouse.Add(closestSegment, go.GetComponent<ArtGalleryLightHouse>());

            CheckSolution();
        }
    }
}