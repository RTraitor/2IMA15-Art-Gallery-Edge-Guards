namespace Util.Geometry {
    using System;
    using System.Collections;
    using System.Collections.Generic;
    using System.Text;
    using UnityEngine;
    using Util.Math;

    /// <summary>
    /// Class representing a line segment consisting of colinear points as well as a begin and end point.
    /// Has auxiliary methods for intersection, orientation and distance among others.
    /// </summary>
    public class MultiLineSegment:LineSegment, IEquatable<MultiLineSegment> {

        // Contains the (Colinear) points of a MultiLineSegment
        // Note that this not an optimal data structure as it has O(n) insertion time (when not inserting at the start/end of the list)
        // Can e.g. be done more efficiently in a Binary Search Tree, but for the small number of points we use this is not needed
        private readonly LinkedList<Vector2> m_points = new LinkedList<Vector2>();

        // Used for quick access
        public LineSegment m_longestSegment { get; private set; }

        public MultiLineSegment(LineSegment s)
            : base(s.Point1, s.Point2)
        {
            m_longestSegment = s;
            m_points.AddFirst(s.Point1);
            m_points.AddFirst(s.Point2);
        }

        public MultiLineSegment(Vector2 a_point1, Vector2 a_point2)
            : this(new LineSegment(a_point1, a_point2)) { }

        public MultiLineSegment(PolarPoint2D a_point1, PolarPoint2D a_point2)
            : this(a_point1.Cartesian, a_point2.Cartesian) { }

        /// <summary>
        /// Adds a point to the multilinesegment.
        /// </summary>
        /// <param name="p">The point to add</param>
        public void AddPoint(Vector2 p) {
            if (MathUtil.EqualsEps(p, Point1) || MathUtil.EqualsEps(p, Point2)) {
                // Point was already in the MultiLineSegment
                return;
            }

            if (!Line.Colinear(p, Point1, Point2)) {
                throw new GeomException("Point " + p + " is not Colinear with the begin ("
                    + Point1 + ") and end (" + Point2 + ") points of the MultiLineSegment");
            }

            // Check if p is between any existing points
            for (LinkedListNode<Vector2> node = m_points.First; node != null; node = node.Next) {
                if (!MathUtil.EqualsEps(node.Value, p)) {
                    if ((new LineSegment(Point1, node.Value)).IsOnSegment(p)) {
                        m_points.AddBefore(node, p);
                        return;
                    }
                } else { // point was already in the point list
                    return;
                }
            }

            // p is not on the current largest line segment
            if ((new LineSegment(Point1, p).IsOnSegment(Point2))) {
                // p is the new end point
                m_points.AddLast(p);
                m_longestSegment = new LineSegment(Point1, p);
            } else {
                // p is the new start point
                m_points.AddFirst(p);
                m_longestSegment = new LineSegment(p, Point2);
            }
        }

        /// <summary>
        /// Allows enumeration over the segment's points.
        /// </summary>
        /// <returns>An enumerable over the segment's points</returns>
        public IEnumerable<Vector2> Points() {
            return m_points;
        }

        /// <summary>
        /// Allows enumeration over the line segments
        /// </summary>
        /// <returns>An enumerable over the line segments</returns>
        public IEnumerable<LineSegment> Segments() {
            List<LineSegment> l = new List<LineSegment>();

            if (m_points.First != null) {
                LinkedListNode<Vector2> v = m_points.First;
                for (LinkedListNode<Vector2> node = m_points.First.Next; node != null; node = node.Next) {
                    l.Add(new LineSegment(v.Value, node.Value));
                    v = node;
                }
            }
            return l;
        }

        public bool Equals(MultiLineSegment other) {
            // Number of points differ
            if (m_points.Count != other.m_points.Count) {
                return false;
            }

            // All points should Epsilon-equal
            using (var e1 = m_points.GetEnumerator())
            using (var e2 = other.m_points.GetEnumerator()) {
                while (e1.MoveNext() && e2.MoveNext()) {
                    if (!MathUtil.EqualsEps(e1.Current, e2.Current)) {
                        // Pair of points did not equal
                        return false;
                    }
                }
            }
            // All points were Epsilon-equal and there were the same number of points
            return true;
        }

        // Auto-generated
        public override int GetHashCode() {
            return -198651050 + EqualityComparer<LinkedList<Vector2>>.Default.GetHashCode(m_points);
        }

        public override string ToString() {
            StringBuilder sb = new StringBuilder("MultiLineSegment: (", m_points.Count);
            foreach (Vector2 p in m_points) {
                sb.Append(p);
                sb.Append(", ");
            }
            sb.Append(")");
            return sb.ToString();
        }
    }
}
