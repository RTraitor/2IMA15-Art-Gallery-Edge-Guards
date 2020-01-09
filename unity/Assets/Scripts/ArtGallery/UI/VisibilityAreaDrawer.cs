namespace ArtGallery {
    using System.Collections.Generic;
    using UnityEngine;
    using UnityEngine.UI;
    using Util.Algorithms.Triangulation;
    using Util.Geometry.DCEL;

    class VisibilityAreaDrawer : MonoBehaviour {
        // line colors for cut lines
        private readonly Color edgeColor = Color.red; //new Color(102f/255, 194f/255, 165f/255);
        private readonly Color vertexColor = Color.blue; //new Color(141f/255, 160f/255, 203f/255);
        private readonly Color faceColor = Color.green; //new Color(252f/255, 141f/255, 98f/255);
        private readonly float VertexRadius = 0.2f;

        // Boolean for whether to display certain parts of the DCEL
        private bool m_displayEdges = false;
        private bool m_displayVertices = false;
        private bool m_displayFaces = false;

        // solution storing the cut lines
        private DCEL m_areas;

        private Transform m_MyTransform;
        private Material m_LineMaterial;

        /// <summary>
        /// Update the solution to be drawn
        /// </summary>
        public DCEL VisibilityAreas {
            get { return m_areas; }
            set { m_areas = value; UpdateDraw(); }
        }

        /// <summary>
        /// Toggles the drawing of the edges of the visiblity areas
        /// </summary>
        public void ToggleDrawEdges() {
            m_displayEdges = !m_displayEdges;
            UpdateDraw();
        }

        /// <summary>
        /// Toggles the drawing of the vertices of the visiblity areas
        /// </summary>
        public void ToggleDrawVertices() {
            m_displayVertices = !m_displayVertices;
            UpdateDraw();
        }

        /// <summary>
        /// Toggles the drawing of the faces of the visiblity areas
        /// </summary>
        public void ToggleDrawFaces() {
            m_displayFaces = !m_displayFaces;
            UpdateDraw();
        }

        /// <summary>
        /// Toggles the drawing of the entire visiblity areas
        /// </summary>
        public void ToggleDrawAll() {
            //m_displayAll = !m_displayAll;
            //UpdateDraw();
            // TODO
        }

        /// <summary>
        /// Updates the line drawing given the divide solution cut lines and toggles.
        /// </summary>
        private void UpdateDraw() {
            //ClearLines();

            //// draw lines that are enabled by adding lines to LineDrawer
            //if (m_displayFaces) {
            //    AddLines(m_solution.Mages, faceColor);
            //}
            //if (m_displayEdges) {
            //    AddLines(m_solution.Archers, edgeColor);
            //}
            //if (m_displayVertices) {
            //    AddLines(m_solution.Spearmen, vertexColor);
            //}
            // TODO
            //if (m_displayAll) {
            //    AddLines(m_solution.All, AllColor);
            //}
        }


        // Use this for initialization
        void Awake() {
            m_MyTransform = this.gameObject.transform;

            // Unity has a built-in shader that is useful for drawing
            // simple colored things.
            Shader shader = Shader.Find("Hidden/Internal-Colored");
            m_LineMaterial = new Material(shader) {
                hideFlags = HideFlags.HideAndDontSave
            };

            // Turn on alpha blending
            m_LineMaterial.SetInt("_SrcBlend", (int) UnityEngine.Rendering.BlendMode.SrcAlpha);
            m_LineMaterial.SetInt("_DstBlend", (int) UnityEngine.Rendering.BlendMode.OneMinusSrcAlpha);

            // Turn backface culling off
            m_LineMaterial.SetInt("_Cull", (int) UnityEngine.Rendering.CullMode.Off);

            // Turn off depth writes
            m_LineMaterial.SetInt("_ZWrite", 0);
        }


        /// <summary>
        /// Draw the edges of the DCEL.
        /// </summary>
        private void DrawEdges() {
            GL.Begin(GL.LINES);

            foreach (HalfEdge edge in m_areas.Edges) {
                GL.Color(edgeColor);
                GL.Vertex3(edge.From.Pos.x, edge.From.Pos.y, 0);
                GL.Vertex3(edge.To.Pos.x, edge.To.Pos.y, 0);
            }
            GL.End();
        }

        /// <summary>
        /// Draws the vertices of the DCEL.
        /// </summary>
        private void DrawVertices() {
            foreach (var vertex in m_areas.Vertices) {
                GL.Begin(GL.TRIANGLE_STRIP);
                GL.Color(vertexColor);

                float step = (2 * Mathf.PI / 200);
                for (float a = 0; a < (2 * Mathf.PI + step); a += step) {
                    //midpoint of the circle.
                    GL.Vertex3(
                        Mathf.Cos(a) * VertexRadius + vertex.Pos.x,
                        Mathf.Sin(a) * VertexRadius + vertex.Pos.y,
                        0f);
                    GL.Vertex(vertex.Pos);
                }
                GL.End();
            }
        }

        /// <summary>
        /// Draws the faces in between the lines.
        /// </summary>
        private void DrawFaces() {
            var triangles = Triangulator.Triangulate(m_areas.InnerFaces).Triangles;

            GL.Begin(GL.TRIANGLES);
            GL.Color(faceColor);

            foreach (var triangle in triangles) {
                GL.Vertex(triangle.P0);
                GL.Vertex(triangle.P1);
                GL.Vertex(triangle.P2);
            }

            GL.End();
        }

        private void OnRenderObject() {
            if (m_areas != null) {
                // Apply the line material
                m_LineMaterial.SetPass(0);

                GL.PushMatrix();
                // Set transformation matrix for drawing to
                // match our transform
                GL.MultMatrix(m_MyTransform.localToWorldMatrix);

                // draw graph components
                DrawFaces();
                DrawEdges();
                DrawVertices();

                GL.PopMatrix();
            }
        }
    }
}
