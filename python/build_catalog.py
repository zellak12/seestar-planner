"""Build targets.csv from OpenNGC (via pyongc). Composite targets use the mean position of their members."""
import sys
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from pyongc import ongc

TARGETS = [  # name, [OpenNGC ids], type
 ("M42_RunningMan", ["NGC1976","NGC1977"], "emission"),
 ("Flame_Horsehead", ["NGC2024","B033"], "emission"),
 ("Pleiades_M45", ["Mel022"], "reflection/cluster"),
 ("California_NGC1499", ["NGC1499"], "emission"),
 ("Rosette_NGC2237", ["NGC2244"], "emission"),
 ("Cone_XmasTree_NGC2264", ["NGC2264"], "emission"),
 ("Seagull_IC2177", ["IC2177"], "emission"),
 ("Jellyfish_IC443", ["IC0443"], "emission (SNR)"),
 ("ThorsHelmet_NGC2359", ["NGC2359"], "emission"),
 ("Beehive_M44", ["NGC2632"], "cluster"),
 ("NGC2903", ["NGC2903"], "galaxy"),
 ("M81_M82", ["NGC3031","NGC3034"], "galaxy"),
 ("LeoTriplet", ["NGC3623","NGC3627","NGC3628"], "galaxy"),
 ("Markarian_Chain", ["NGC4374","NGC4406","NGC4435","NGC4438","NGC4461","NGC4473","NGC4477"], "galaxy"),
 ("M106", ["NGC4258"], "galaxy"),
 ("Needle_NGC4565", ["NGC4565"], "galaxy"),
 ("Whale_Hockey_NGC4631", ["NGC4631","NGC4656"], "galaxy"),
 ("Whirlpool_M51", ["NGC5194"], "galaxy"),
 ("Pinwheel_M101", ["NGC5457"], "galaxy"),
 ("CentaurusA_NGC5128", ["NGC5128"], "galaxy"),
 ("M13", ["NGC6205"], "cluster"),
 ("Antares_M4", ["NGC6121"], "cluster/reflection"),
 ("CatsPaw_NGC6334", ["NGC6334"], "emission"),
 ("Lagoon_Trifid", ["NGC6523","NGC6514"], "emission"),
 ("StarCloud_M24", ["IC4715"], "star cloud"),
 ("Eagle_M16", ["NGC6611"], "emission"),
 ("Swan_M17", ["NGC6618"], "emission"),
 ("Dumbbell_M27", ["NGC6853"], "planetary"),
 ("Crescent_NGC6888", ["NGC6888"], "emission"),
 ("Sadr_IC1318", ["IC1318"], "emission"),
 ("Veil_center", ["NGC6960","NGC6992"], "emission (SNR)"),
 ("NorthAmerica_NGC7000", ["NGC7000"], "emission"),
 ("ElephantsTrunk_IC1396", ["IC1396"], "emission"),
 ("Cocoon_IC5146", ["IC5146"], "emission"),
 ("Helix_NGC7293", ["NGC7293"], "planetary"),
 ("Bubble_M52", ["NGC7635","NGC7654"], "emission"),
 ("NGC55", ["NGC0055"], "galaxy"),
 ("M31", ["NGC0224"], "galaxy"),
 ("NGC253", ["NGC0253"], "galaxy"),
 ("Pacman_NGC281", ["NGC0281"], "emission"),
 ("M33", ["NGC0598"], "galaxy"),
 ("DoubleCluster", ["NGC0869","NGC0884"], "cluster"),
 ("Heart_IC1805", ["IC1805"], "emission"),
 ("Soul_IC1848", ["IC1848"], "emission"),
 ("NGC1097", ["NGC1097"], "galaxy"),
 ("Fornax_NGC1399", ["NGC1399"], "galaxy"),
 ("NGC1365", ["NGC1365"], "galaxy"),
 ("FlamingStar_Tadpoles", ["IC0405","IC0410"], "emission"),
 ("M3", ["NGC5272"], "cluster"),
 ("Sombrero_M104", ["NGC4594"], "galaxy"),
 ("SouthernPinwheel_M83", ["NGC5236"], "galaxy"),
 ("BlackEye_M64", ["NGC4826"], "galaxy"),
 ("Sunflower_M63", ["NGC5055"], "galaxy"),
 ("Antennae_NGC4038", ["NGC4038","NGC4039"], "galaxy"),
 ("M5", ["NGC5904"], "cluster"),
 ("Ring_M57", ["NGC6720"], "planetary"),
 ("WildDuck_M11", ["NGC6705"], "cluster"),
 ("M92", ["NGC6341"], "cluster"),
 ("Butterfly_M6", ["NGC6405"], "cluster"),
 ("Ptolemy_M7", ["NGC6475"], "cluster"),
 ("Prawn_IC4628", ["IC4628"], "emission"),
 ("Fireworks_NGC6946", ["NGC6946"], "galaxy"),
 ("M71", ["NGC6838"], "cluster"),
]

out = open(sys.argv[1], "w")
out.write("# J2000 positions from the OpenNGC catalog (github.com/mattiaverga/OpenNGC) via pyongc.\n")
out.write("# Composite targets (pairs/groups) use the mean position of the listed members.\n")
out.write("# Add your own lines: name,ra,dec[,type,members]  (ra hh:mm:ss, dec ±dd:mm:ss; or decimal degrees)\n")
out.write("name,ra,dec,type,members\n")
for name, ids, typ in TARGETS:
    cs = []
    for i in ids:
        o = ongc.get(i)
        if o is None or o.ra is None:
            raise SystemExit(f"{name}: {i} not found / no coords")
        cs.append(SkyCoord(o.ra, o.dec, unit=(u.hourangle, u.deg)))
    v = np.mean([c.cartesian.xyz.value for c in cs], axis=0)
    m = SkyCoord(x=v[0], y=v[1], z=v[2], representation_type="cartesian").spherical
    c = SkyCoord(m.lon, m.lat)
    ra = c.ra.to_string(u.hourangle, sep=":", precision=0, pad=True)
    dec = c.dec.to_string(u.deg, sep=":", precision=0, alwayssign=True, pad=True)
    out.write(f"{name},{ra},{dec},{typ},{' '.join(ids)}\n")
print("ok")
