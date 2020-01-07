namespace ArtGallery {

    using System.Collections.Generic;
    using UnityEngine;
    using Util.Geometry;
    using Util.Geometry.DCEL;
    using Util.Geometry.Polygon;
    using Main;
    using System.Linq;

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
                    DCELVertex vertex = ((List<DCELVertex>) faceID.Value.OuterVertices)[0];
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
        /// Computes the weakly visible points in polygon P from point z
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
                    if (u0 == null)
                    {
                        u0 = intersection;
                    } else
                    {
                        if (intersection?.x < u0?.x)
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
            Polygon2D poly = new Polygon2D(vertices); //ADDED, NOT SURE IF NECESSARY

            // Initially the stack contains u0 and u1
            visiblePoints.Push((Vector2) u0);
            visiblePoints.Push(vertices[1]);
            // Next element to check is u2
            int nextElement = vertices.IndexOf(vertices[2]);
            // Algorithm terminates when u0 is pushed on the stack again
            while (!nextElement.Equals(vertices.IndexOf((Vector2) u0)))
            {
                // Get index of element u_{i-2}
                int element2Prev = nextElement - 2; 
                
                //Check if u_{i-2} lies to the right or left of line segment (z, s_j)
                if (new LineSegment(z, visiblePoints.Peek()).IsRightOf(vertices[element2Prev]))
                {
                    // Case C1
                } 
                else
                {
                    // Case C2
                }
            }

            return visiblePoints.ToList();
        }

        /// <summary>
        /// Approximate the minimum number of edge guards needed for the current level's polygon
        /// with approximation ratio O(log(n))
        /// </summary>
        /// <returns>
        /// An approximation of the minimum number of edge guards needed for the current level's 
        /// polygon with approximation ratio O(log(n))
        /// </returns>
        private int calcEdgeGuards()
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
            // return if lighthouse was already selected or player can place no more lighthouses
            if (m_selectedLighthouse != null || m_solution.Count >= m_maxNumberOfLighthouses)
                return;

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