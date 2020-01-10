using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
namespace Main
{
    public class SetEqualityComparer : IEqualityComparer<HashSet<int>>
    {
        public bool Equals(HashSet<int> x, HashSet<int> y) => x.SetEquals(y);

        public int GetHashCode(HashSet<int> obj)
        {
            var code = 0;
            foreach (var value in obj)
            {
                code += value;
            }
            return code;
        }
    }

    public class SetSizeComparer : IComparer<KeyValuePair<int, HashSet<int>>>
    {
        public int Compare(KeyValuePair<int, HashSet<int>> x, KeyValuePair<int, HashSet<int>> y)
        {
            bool testLess = x.Value.Count() < y.Value.Count();
            bool testEqual = x.Value.Count() == y.Value.Count();
            if (testLess)
            {
                // descending order
                return 1;
            }
            else if (testEqual)
            {
                return 0;
            }
            return -1;
        }
    }

    public class SetCover
    {
        /// <summary>
        /// Remove all duplicate hashsets
        /// </summary>
        /// <param name="S"></param>
        /// <returns>A dictionary without duplicates</returns>
        private static Tuple<bool, Dictionary<int, HashSet<int>>> RemoveDuplicates(Dictionary<int, HashSet<int>> S)
        {
            var distinct = S.Select(x => x.Value).Distinct(new SetEqualityComparer());
            if (distinct.Count() == S.Count())
            {
                return Tuple.Create(false, S);
            }
            else
            {
                var result = new Dictionary<int, HashSet<int>>();
                foreach (var set in distinct)
                {
                    for (int i = 0; i < S.Count(); i++)
                    {
                        var pair = S.ElementAt(i);
                        if (pair.Value.SetEquals(set))
                        {
                            result.Add(pair.Key, set);
                            break;
                        }
                    }
                }
                return Tuple.Create(true, result);
            }
        }

        /// <summary>
        /// Reduce an instance of Set Cover using safe reduction rules
        /// </summary>
        /// <param name="U">The universe</param>
        /// <param name="F">The familiy of sets</param>
        /// <returns>
        /// (U, F, delta), where U and F are "irreducible" and
        /// delta is the number of sets that must be included
        /// </returns>
        private static Tuple<HashSet<int>, Dictionary<int, HashSet<int>>, List<int>> Reduce(HashSet<int> U, Dictionary<int, HashSet<int>> F)
        {
            var FaceSets = new Dictionary<int, HashSet<int>>();
            foreach (var key in F.Keys)
            {

                foreach (var face in F[key])
                {
                    if (FaceSets.ContainsKey(face))
                    {
                        FaceSets[face].Add(key);
                    }
                    else
                    {
                        FaceSets.Add(face, new HashSet<int> { key });
                    }
                }
            }


            // main reduction loop
            // whether the solution changed at all
            var reduced = false;

            // whether the solution changed this iteration
            var changed = true;
            var delta = new List<int>();
            do
            {
                changed = false;
                // last face was removed 
                if (FaceSets.Count() == 0)
                {
                    return Tuple.Create(new HashSet<int>(), new Dictionary<int, HashSet<int>>(), delta);
                }

                // remove duplicates
                var noDups = RemoveDuplicates(F);
                if (noDups.Item1)
                {
                    reduced = true;
                    changed = true;
                    F = noDups.Item2;
                    continue;
                }

                // remove isolated faces
                var isolated =
                from set in FaceSets
                where (set.Value.Count == 1)
                select set.Key;

                if (isolated.Count() > 0)
                {
                    reduced = true;
                    changed = true;

                    // collect all faces of sets which include the isolated faces                    
                    var coveredFaces = new HashSet<int>();
                    var coveredEdges = new HashSet<int>();
                    foreach (var face in isolated)
                    {
                        delta.Add(face);
                        var ids = FaceSets[face];
                        // only executed once...
                        foreach (var id in ids)
                        {
                            foreach (var f in F[id])
                            {
                                coveredFaces.Add(f);
                            }
                            coveredEdges.Add(id);
                        }
                    }

                    // build new FaceSets
                    var newFaceSets = new Dictionary<int, HashSet<int>>();
                    foreach (var pair in FaceSets)
                    {
                        if (!coveredFaces.Contains(pair.Key))
                        {
                            var newSet = pair.Value;
                            newSet.ExceptWith(coveredEdges);
                            newFaceSets.Add(pair.Key, newSet);
                        }
                    }
                    FaceSets = newFaceSets;

                    // build new F
                    var newF = new Dictionary<int, HashSet<int>>();
                    foreach (var pair in F)
                    {
                        if (!coveredEdges.Contains(pair.Key))
                        {
                            var newSet = pair.Value;
                            newSet.ExceptWith(coveredFaces);
                            newF.Add(pair.Key, newSet);
                        }
                    }
                    F = newF;

                    continue;
                }


                // find the first completely included subsets and remove it
                for (int i = 0; i < F.Count && !changed; i++)
                {
                    var key = F.ElementAt(i).Key;
                    var A = F.ElementAt(i).Value;
                    for (int j = i + 1; j < F.Count; j++)
                    {
                        bool complete = true;
                        var B = F.ElementAt(j).Value;
                        foreach (var x in A)
                        {
                            if (!B.Contains(x))
                            {
                                complete = false;
                            }
                        }
                        if (complete)
                        {
                            foreach (var x in A)
                            {
                                FaceSets[x].Remove(key);
                            }
                            F.Remove(key);

                            reduced = true;
                            changed = true;
                            break;
                        }

                    }
                }

                if (changed)
                {
                    continue;
                }

                // remove the first face which is present in every set
                var included =
                from set in FaceSets
                where (set.Value.Count == F.Count)
                select set.Key;

                if (included.Count() > 0)
                {
                    reduced = true;
                    changed = true;
                    var face = included.ElementAt(0);

                    FaceSets.Remove(face);

                    foreach (var set in F)
                    {
                        set.Value.Remove(face);
                    }
                    continue;
                }

            } while (changed);


            // recover the instance
            if (reduced)
            {

                var new_U = new HashSet<int>();
                foreach (var face in FaceSets)
                {
                    new_U.Add(face.Key);
                }

                return Tuple.Create(new_U, F, delta);
            }
            else
            {
                return Tuple.Create(U, F, delta);
            }
        }

        /// <summary>
        /// Generate all possible tuples (of indices) of size k
        /// </summary>
        /// <param name="S">list of all generated tuples</param>
        /// <param name="A">current list of indices part of a tuple</param>
        /// <param name="j">index so far</param>
        /// <param name="k">size of the tuples</param>
        /// <param name="n">size of A</param>
        /// <returns>S</returns>
        private static List<int[]> GenerateTuples(List<int[]> S, int[] A, int j, int k, int n)
        {

            if (k != 0)
            {
                // first index has no predecessor
                if (j == 0)
                {
                    for (var i = 0; i < n; i++)
                    {
                        A[0] = i;
                        GenerateTuples(S, A, j + 1, k - 1, n);
                    }
                }
                else
                {
                    // pick the next index compared to its predecessor
                    for (var i = A[j - 1] + 1; i < n; i++)
                    {
                        A[j] = i;
                        GenerateTuples(S, A, j + 1, k - 1, n);
                    }
                }
            }
            else
            {
                S.Add(A);
            }
            // note that S is not copied at each recursive call
            // instead, S is a reference that can be updated
            return S;
        }

        /// <summary>
        /// Find the exact (optimal) solution of the instance (U, F)
        /// </summary>
        /// <param name="U">The universe</param>
        /// <param name="F">The familiy of subsets</param>
        /// <returns>The optimal solution of the instance (U, F)</returns>
        private static List<int> ExactSetCover(HashSet<int> U, Dictionary<int, HashSet<int>> F)
        {
            var result = new List<int>();
            // Linear search, could be improved to binary search
            for (int k = 1; k <= F.Count(); k++)
            {
                // Initialization for tuple generation
                var tuples = new List<int[]>();
                var A = new int[F.Count()];
                var j = 0;

                // Get all k-tuples
                tuples = GenerateTuples(tuples, A, j, k, F.Count());
                foreach (var tuple in tuples)
                {
                    // set representation
                    var S = new HashSet<int>();
                    // compute union of selected sets
                    for (var i = 0; i < k; i++)
                    {
                        foreach (var item in F.ElementAt(tuple[i]).Value)
                        {
                            S.Add(item);
                        }
                    }

                    // Check if the sets are equal
                    if (S.SetEquals(U))
                    {
                        for (int i = 0; i < k; i++)
                        {
                            result.Add(tuple[i]);
                        }
                        return result;
                    }
                }
            }
            // unreachable as there is always a solution, namely all edges.
            return result;
        }

        /// <summary>
        /// Greedily solve Set Cover
        /// </summary>
        /// <param name="U">The universe</param>
        /// <param name="F">The family of subsets</param>
        /// <returns>The selected sets by id that are a solution</returns>
        private static List<int> GreedySetCover(HashSet<int> U, Dictionary<int, HashSet<int>> F)
        {
            var result = new List<int>();
            var listF = F.ToList();


            while (U.Count != 0)
            {
                var maxSize = 0;
                var max = new KeyValuePair<int, HashSet<int>>();
                foreach (var pair in F)
                {
                    if (pair.Value.Count() > maxSize)
                    {
                        maxSize = pair.Value.Count;
                        max = pair;
                    }
                }

                result.Add(max.Key);
                
                // safe set subtraction
                U.ExceptWith(max.Value);
                foreach (var pair in F)
                {
                    pair.Value.ExceptWith(max.Value);
                }
            }
            return result;
        }

        /// <summary>
        /// Method to be called to solve an instance (U, F)
        /// </summary>
        /// <param name="U">The universe</param>
        /// <param name="F">Yhe familiy of subsets</param>
        /// <returns>
        /// The size of a solution of the instance (U, F), optimal or 
        /// estimate.
        /// </returns>
        public static List<int> Solve(HashSet<int> U, Dictionary<int, HashSet<int>> F)
        {
            // check if input is within bounds
            if (F.Count <= 25 && U.Count <= 500)
            {
                // reduce it
                var t = Reduce(U, F);
                U = t.Item1;
                F = t.Item2;
                var delta = t.Item3;

                if (U.Count == 0)
                {
                    return delta;
                }
                else
                {
                    var result = ExactSetCover(U, F);
                    foreach (var id in result)
                    {
                        delta.Add(id);
                    }
                    return delta;
                }
            }
            else
            {
                return GreedySetCover(U, F);
            }
        }
    }
}
