
namespace ArtGallery {
    using System.Collections.Generic;
    using Util.Geometry;
    using Util.Math;

    internal class UndirectedSegmentComparer : IEqualityComparer<LineSegment> {
        public UndirectedSegmentComparer() {}

        public bool Equals(LineSegment x, LineSegment y) {
            return (MathUtil.EqualsEpsVertex(x.Point1, y.Point1) &&
                        MathUtil.EqualsEpsVertex(x.Point2, y.Point2))
                    || (MathUtil.EqualsEpsVertex(x.Point1, y.Point2) &&
                        MathUtil.EqualsEpsVertex(x.Point2, y.Point1));
        }

        public int GetHashCode(LineSegment l) {
            return l.GetHashCode();
        }
    }
}