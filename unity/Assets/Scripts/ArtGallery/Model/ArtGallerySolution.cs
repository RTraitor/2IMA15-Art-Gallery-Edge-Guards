namespace ArtGallery
{
    using System.Collections.Generic;
    using System.Linq;
    using UnityEngine;
    using Util.Algorithms.Polygon;
    using Util.Geometry.Polygon;

    /// <summary>
    /// Stores a list of lighthouses and can calculate their combined visibility area.
    /// Handles destruction of lighthouse game objects.
    /// </summary>
    public class ArtGallerySolution : ScriptableObject
    {

        /// <summary>
        /// The number of lighthouses placed
        /// </summary>
        public int Count {  get { return m_lighthouses.Count; } }

        /// <summary>
        /// Total area visible by all lighthouses
        /// </summary>
        public float Area
        {
            get
            {
                if (Count <= 0) return 0f;

                // create multi polygon of visibility area
                var visiblePolygon = new MultiPolygon2D(m_lighthouses[0].VisionPoly);

                // add visibility polygons, cutting out the overlap
                foreach (ArtGalleryLightHouse lighthouse in m_lighthouses.Skip(1))
                {
                    visiblePolygon = Clipper.CutOut(visiblePolygon, lighthouse.VisionPoly);
                    visiblePolygon.AddPolygon(lighthouse.VisionPoly);
                }

                // add visibility polygons for edge guards
                foreach (ArtGalleryLightHouseInvis invis in m_lighthousesinvis)
                {
                    visiblePolygon = Clipper.CutOut(visiblePolygon, invis.VisionPoly);
                    visiblePolygon.AddPolygon(invis.VisionPoly);
                }

                // return total area
                return visiblePolygon.Area;
            }
        }

        // collection of lighthouses
        private List<ArtGalleryLightHouse> m_lighthouses;

        // collection of invis lighthouses (EDGE GUARDS ONLY)
        private List<ArtGalleryLightHouseInvis> m_lighthousesinvis;

        // stores lighthouse objects for easy destroyal
        private List<GameObject> m_objects;

        public ArtGallerySolution()
        {
            m_objects = new List<GameObject>();
            m_lighthouses = new List<ArtGalleryLightHouse>();
            m_lighthousesinvis = new List<ArtGalleryLightHouseInvis>();
        }

        /// <summary>
        /// Add lighthouse to solution.
        /// </summary>
        /// <param name="m_lighthouse"></param>
        public void AddLighthouse(ArtGalleryLightHouse m_lighthouse)
        {
            m_lighthouses.Add(m_lighthouse);
        }

        /// <summary>
        /// Add lighthouse invis to solution. (EDGE GUARDS ONLY)
        /// </summary>
        /// <param name="m_invis"></param>
        public void AddLighthouseInvis(ArtGalleryLightHouseInvis m_invis)
        {
            m_lighthousesinvis.Add(m_invis);
        }

        /// <summary>
        /// Create a lighthouse object for given game object and add to solution.
        /// </summary>
        /// <param name="obj"></param>
        public void AddLighthouse(GameObject obj)
        {
            // remember object for removal
            m_objects.Add(obj);

            // Check if we add invis lighthouses or normal lighthouses ( relevant for edge guards)
            if (obj.GetComponent<ArtGalleryLightHouse>() == null)
            {
                AddLighthouseInvis(obj.GetComponent<ArtGalleryLightHouseInvis>());
            }
            else
            {
                // add the lighthouse component of game object to solution
                AddLighthouse(obj.GetComponent<ArtGalleryLightHouse>());

            }
        }

        /// <summary>
        /// Remove the given lighthouse from the solution
        /// </summary>
        /// <param name="m_lighthouse"></param>
        public void RemoveLighthouse(ArtGalleryLightHouse m_lighthouse)
        {
            m_lighthouses.Remove(m_lighthouse);
            Destroy(m_lighthouse);
        }

        /// <summary>
        /// Remove the given lighthouse invis from the solution (EDGE GUARDS ONLY)
        /// </summary>
        /// <param name="m_lighthouse"></param>
        public void RemoveLighthouseInvis(ArtGalleryLightHouseInvis m_invis)
        {
            m_lighthousesinvis.Remove(m_invis);
            Destroy(m_invis);
        }

        /// <summary>
        /// Clears the lighthouse lists and destroys all corresponding game objects
        /// </summary>
        public void Clear()
        {
            foreach (var lh in m_lighthouses) Destroy(lh);
            //for edge guards
            foreach (var invis in m_lighthousesinvis) Destroy(invis);
            foreach (var obj in m_objects) Destroy(obj);
            m_lighthouses.Clear();
            //For edge guards
            m_lighthousesinvis.Clear();
        }

        public void OnDestroy()
        {
            Clear();
        }
    }
}
