namespace ArtGallery
{
    using UnityEngine;
    using Util.Geometry.Polygon;
    using Util.Geometry;

    /// <summary>
    /// Represents invisible lighthouse objects in the game,
    /// to save the visibility polygons with edge guards.
    /// </summary>
    public class ArtGalleryLightHouseInvis : MonoBehaviour
    {
        // stores a prefab object for the vision polygon
        [SerializeField]
        private GameObject m_visionAreaPrefab;

        private ArtGalleryController m_controller;

        public LineSegment m_segment;

        /// <summary>
        /// Mesh variable of the art gallery.
        /// </summary>
        public ArtGalleryIsland VisionAreaMesh { get; set; }

        /// <summary>
        /// Stores lighthouse position. Updates vision after a change in position.
        /// </summary>
        /*
        public Vector3 Pos
        {
            get
            {
                return gameObject.transform.position;
            }
            set
            {
                gameObject.transform.position = value;

                // update vision polygon
                m_controller.UpdateVision(this);
            }
        }
        */

        /// <summary>
        /// Holds the visibility polygon.
        /// </summary>
        public Polygon2D VisionPoly { get; set; }

        // Use this for initialization
        void Awake()
        {
            m_controller = FindObjectOfType<ArtGalleryController>();

            // initialize the vision polygon
            GameObject go = Instantiate(m_visionAreaPrefab, new Vector3(0, 0, -1.5f), Quaternion.identity) as GameObject;
            VisionAreaMesh = go.GetComponent<ArtGalleryIsland>();

            //UpdateVision();
        }

        void OnDestroy()
        {
            if (VisionAreaMesh != null)
            {
                Destroy(VisionAreaMesh.gameObject);
            }
        }

        /*
        void OnMouseDown()
        {
            m_controller.SelectLighthouse(this);
        }

        public void UpdateVision()
        {
            m_controller.UpdateVision(this);
        }
        */
    }
}
