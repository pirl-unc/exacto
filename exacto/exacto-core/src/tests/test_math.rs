// use crate::common::math::{betabinom_sf, betabinom_sf_complement};


// /// `betabinom_sf_complement` and `betabinom_sf` are the same function computed from
// /// opposite tails; they must agree wherever both are well-conditioned.
// #[test]
// fn test_betabinom_sf_complement_matches_direct_sum() {
//     for &n in &[1u64, 2, 10, 100, 1_000] {
//         for &k in &[0u64, 1, 2, 3, 5, 10] {
//             if k > n {
//                 continue;
//             }
//             for &(eps, rho) in &[(0.01, 0.01), (0.001, 0.02), (0.1, 0.1)] {
//                 let direct: f64 = betabinom_sf(k, n, eps, rho);
//                 let complement: f64 = betabinom_sf_complement(k, n, eps, rho);
//                 assert!(
//                     (direct - complement).abs() <= 1e-9,
//                     "sf({k}, {n}, {eps}, {rho}): direct {direct} vs complement {complement}"
//                 );
//             }
//         }
//     }
// }
