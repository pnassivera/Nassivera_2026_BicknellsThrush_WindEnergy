# Supplementary methodological note: rasterization of infrastructure footprints

Habitat loss and fragmentation were assessed using the R package `landscapemetrics`, which calculates landscape metrics from categorical rasters. Each cell is assigned to a single habitat class, requiring a decision on how to represent cells partially intersected by infrastructure.

The species distribution model (SDM) used in this study has a **10-m spatial resolution**, providing a fine-scale representation of habitat distribution across the landscape.

To accommodate the SDM's spatial uncertainty, a **10-m buffer** was applied around infrastructure footprints before rasterization. The buffered footprints were then rasterized using `terra::rasterize()` with **`touches = TRUE`**, meaning that any cell intersected by the buffered footprint was classified as affected, regardless of the proportion of its area covered.

This approach was chosen to avoid overlooking potentially affected habitat cells and the disruption of narrow connections between habitat patches, particularly where linear infrastructure crosses habitat corridors. It is also consistent with a precautionary conservation perspective, as ecological disturbance may extend beyond the physical infrastructure footprint due to edge effects (Pfeifer *et al.*, 2017).

**Consequently, habitat loss estimates reported in Nassivera *et al.* (2026) are inherently limited by the 10-m spatial resolution of the SDM, and should be interpreted as estimates of potentially affected habitat at the resolution of the SDM, rather than exact measurements of physical habitat removal or ecological disturbance.**
