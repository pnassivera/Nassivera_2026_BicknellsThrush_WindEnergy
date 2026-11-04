# Supplementary methodological note: rasterization of infrastructure footprints

Habitat loss and fragmentation were assessed using the R package `landscapemetrics`, which calculates landscape metrics from categorical rasters. Each cell is assigned to a single habitat class, requiring a decision on how to represent cells partially intersected by infrastructure.

The species distribution model (SDM) used in this study has a **10-m spatial resolution**, providing a fine-scale representation of habitat distribution across the landscape.

To accommodate the SDM's spatial uncertainty, a **10-m geometric buffer** was applied around infrastructure footprints before rasterization. The buffered footprints were then rasterized using `terra::rasterize()` with **`touches = TRUE`**, meaning that any cell intersected by the buffered footprint was classified as affected, regardless of the proportion of its area covered.

This conservative approach was chosen to avoid overlooking potentially affected habitat cells and the disruption of narrow connections between habitat patches, particularly where linear infrastructure crosses habitat corridors. It is also consistent with a precautionary conservation perspective, as ecological disturbance may extend beyond the physical infrastructure footprint.

**Consequently, the reported results should be interpreted as conservative estimates of affected habitat at the 10-m resolution of the SDM, rather than exact measurements of physical habitat removal or ecological disturbance.**
