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

        //Unity references
        private VisibilityAreaDrawer m_areaDrawer;

        // Update is called once per frame COULD GIVE PROBLEMS - TEST NECESSARY
        void Update() {
            //handle input key presses
            if (Input.GetKeyDown("a")) {
                m_areaDrawer.ToggleDrawAll();
            }
            if (Input.GetKeyDown("s")) {
                m_areaDrawer.ToggleDrawEdges();
            }
            if (Input.GetKeyDown("v")) {
                m_areaDrawer.ToggleDrawVertices();
            }
            if (Input.GetKeyDown("f")) {
                m_areaDrawer.ToggleDrawFaces();
            }
        }

        public override void InitLevel() {
            base.InitLevel();
            Debug.Log("Initialising Level...");
            Debug.Log(LevelPolygon.Segments.ToString());
            DCEL dcell = new DCEL();
            foreach (LineSegment s in LevelPolygon.Segments) {
                Debug.Log("Segment: " + s.ToString());
                dcell.AddSegment(s);
            }

            m_areaDrawer.VisibilityAreas = dcell;
        }

        /// <summary>
        /// Computes the visibilty regions of a given polygon and outputs it as a DCEL
        /// </summary>
        /// <param name="poly">The polygon to compute the visibility regions for</param>
        /// <returns>
        /// A DCEL containing the visibility regions of a polygon
        /// </returns>
        private DCEL computeVisibilityRegions(Polygon2D poly) {
            ICollection<LineSegment> segments = getVisibilitySegments(poly);
            return mergeSegments(segments);
        }

        /// <summary>
        /// Gets the (overlapping) Line Segments that are needed to compute a given polygon's visibility regions
        /// </summary>
        /// <param name="poly">The polygon to compute the segments for</param>
        /// <returns>
        /// The (overlapping) Line Segments that are needed to compute the given polygon's visbility regions
        /// </returns>
        private ICollection<LineSegment> getVisibilitySegments(Polygon2D poly) {
            LinkedList<LineSegment> segments = new LinkedList<LineSegment>();
            
            foreach (Vector2 v in poly.Vertices) {
                foreach (Vector2 v2 in poly.Vertices) {
                    if (v.Equals(v2)) {
                        continue;
                    }

                    if (poly.isConvex(v2) == false) { // v2 is Reflex

                        LineSegment s = new LineSegment(v, v2);
                        if (isIntersectPolygon(poly, s)) {
                            // Skip if the two vertices cannot see each other.
                            // NOTE: This includes colinear points (meaning that only 'adjacent' colinear points are not skipped)
                            continue;
                        }

                        Vector2? closestIntersection = intersectPolygonClosest(poly, new Line(v, v2));

                        if (closestIntersection != null) {
                            segments.AddFirst(new LineSegment(v2, (Vector2) closestIntersection));
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
        private Vector2? intersectPolygonClosest(Polygon2D poly, Line l) {
            LineSegment smallestSegment = null;
            foreach (LineSegment s in poly.Segments) {
                // Ignore the segments that the ray starts in (if any)
                if (s.IsOnSegment(l.Point2)) {
                    continue;
                }

                // Will contain the intersection of l and s (if any)
                Vector2? intersection = null;

                // Edge Case: Handle segments that completely overlap the given line
                if (l.IsOnLine(s.Point1) && l.IsOnLine(s.Point2)) {
                    // All points of the segment overlap!
                    if ((new LineSegment(l.Point2, s.Point1)).IsOnSegment(s.Point2)) {
                        // Endpoint of the segment is closer
                        intersection = s.Point2;
                    } else if ((new LineSegment(l.Point2, s.Point2)).IsOnSegment(s.Point1)) {
                        // Beginpoint of the segment is closer
                        intersection = s.Point1;
                    }
                } else {
                    // A single point of the segment can overlap
                    intersection = s.Intersect(l);
                }

                // Ignore intersections on the 'wrong' side
                // TODO: optimize this(?); ignore segments on one side of the line perpendicular to l
                if (intersection != null && (new LineSegment(l.Point1, (Vector2) intersection)).IsOnSegment(l.Point2)) {
                    continue;
                }

                // There is an intersection, and it was closer than the previous closest intersection
                if (intersection != null && (smallestSegment == null || smallestSegment.IsOnSegment((Vector2) intersection))) {
                    smallestSegment = new LineSegment(l.Point2, (Vector2) intersection);
                }
            }
            if (smallestSegment != null) {
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
        private bool isIntersectPolygon(Polygon2D poly, LineSegment vertexSegment) {
            foreach (LineSegment s in poly.Segments) {
                // Skip the segments that have a begin or end point in common with the given segment
                if (s.Point1 == vertexSegment.Point1 || s.Point2 == vertexSegment.Point2
                        || s.Point1 == vertexSegment.Point2 || s.Point2 == vertexSegment.Point1) {
                    continue;
                }

                // Edge Case: Skip segments whose begin and end points are BOTH colinear with the given segment
                //if (Line.Colinear(s.Point1, s.Point2, vertexSegment.Point1)
                //        && Line.Colinear(s.Point1, s.Point2, vertexSegment.Point2)) {
                //    continue;
                //}

                // Skip begin and end points of segments
                Vector2? x = vertexSegment.Intersect(s);
                if (x != null) {
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
        private DCEL mergeSegments(ICollection<LineSegment> segments) {
            
            List<LineSegment> lineSegments = new List<LineSegment>();
            // For each pair of segments
            foreach (LineSegment s1 in segments) {
                foreach (LineSegment s2 in segments) {
                    // Segments must not be the same
                    if (!s1.Equals(s2)) {
                        Vector2? intersection = s1.Intersect(s2);
                        
                        // TODO: Handle multiple intersections in the same linesegment

                        // If the segments intersect, create a new vertex at the intersecion and split the segments
                        if (intersection != null) {
                            lineSegments.Add(new LineSegment(s1.Point1, (Vector2) intersection));
                            lineSegments.Add(new LineSegment(s1.Point2, (Vector2) intersection));
                            lineSegments.Add(new LineSegment(s2.Point1, (Vector2) intersection));
                            lineSegments.Add(new LineSegment(s2.Point2, (Vector2) intersection));
                        }
                    }
                }
            }

            // Create the DCEL
            DCEL dcel = new DCEL();

            // Add the segments to the DCEL
            foreach (LineSegment segment in lineSegments) {
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
            //// Draw lines through every pair of vertices 
            //List<Line> lines = generateLines((List<Vector2>) LevelPolygon.Vertices);
            //// Create convex components 
            //List<LineSegment> lineSegments = generateLineSegments(lines);
            //DCEL dcel = createDCELFromPolygonAndSegments(LevelPolygon, lineSegments);
            //List<Face> faces = (List<Face>) dcel.InnerFaces;
            //// Give IDs to the faces
            //setFaceIDs(faces);
            //// Give IDs to the edges
            //setEdgeIDs((List<LineSegment>) LevelPolygon.Segments);
            //// Store visible convex components(/faces) for each edge by index
            //computeVisibleComponentsPerEdge(dcel);
            //// Calculate the needed number of edge guards
            //return SetCover.Solve(new HashSet<int>(faceIDs.Keys), visibleCompIDsPerEdgeID);
            return 5;
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