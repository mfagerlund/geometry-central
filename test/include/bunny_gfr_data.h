#pragma once

// C# GFR comparison data for Bunny mesh
// From Bunny_comparison.json - IterativeDistance values

struct BunnyGFRPair {
  size_t from;
  size_t to;
  double csharpGfrDistance;  // C# GFR (IterativeDistance)
  double flipOutDistance;     // FlipOut reference
};

// First 50 pairs for quick validation
static const BunnyGFRPair BUNNY_GFR_PAIRS[] = {
    {980, 1526, 3.2585341930389404, 3.25941985083333},
    {1238, 2481, 1.4036065340042114, 1.41471780673263},
    {908, 1405, 1.833968162536621, 1.83570895406115},
    {1987, 1102, 1.0388470888137817, 1.03884708274685},
    {982, 505, 1.7909395694732666, 1.7909395900885},
    {577, 1643, 3.876040458679199, 3.87604043192852},
    {1656, 640, 2.268129587173462, 2.26812970322917},
    {385, 1576, 0.33970382809638977, 0.339703832317581},
    {469, 2449, 2.3622918128967285, 2.36520716456929},
    {419, 1679, 3.030559778213501, 3.03055988582219},
};
static const size_t BUNNY_GFR_COUNT = 10;
