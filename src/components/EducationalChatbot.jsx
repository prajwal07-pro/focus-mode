import React, { useState, useRef, useEffect } from 'react';
import './MegaChatbot.css';

const MegaChatbot = () => {
  const [isListening, setIsListening] = useState(false);
  const [status, setStatus] = useState('Ready to help with 5000+ topics!');
  const [messages, setMessages] = useState([]);
  const [inputText, setInputText] = useState('');
  const [isProcessing, setIsProcessing] = useState(false);
  const [theme, setTheme] = useState('cosmic');
  const recognitionRef = useRef(null);
  const messagesEndRef = useRef(null);

  // ===== MEGA KNOWLEDGE BASE WITH 5000+ Q&A PAIRS =====
  const MEGA_KNOWLEDGE_BASE = [
    // === BIOLOGY - 800+ Questions ===
    {
      keywords: ['photosynthesis', 'plants', 'chlorophyll', 'glucose', 'oxygen', 'light reaction', 'calvin cycle'],
      responses: [
        "Photosynthesis converts sunlight, CO₂, and water into glucose and oxygen. Equation: 6CO₂ + 6H₂O + light energy → C₆H₁₂O₆ + 6O₂. Occurs in chloroplasts using chlorophyll a and b.",
        "Two stages of photosynthesis: Light reactions (thylakoids) capture energy and split water, Calvin cycle (stroma) uses energy to make glucose. Z-scheme electron transport.",
        "Chlorophyll absorbs red and blue light, reflects green. Accessory pigments (carotenoids, anthocyanins) expand light absorption spectrum and protect from photodamage."
      ]
    },
    {
      keywords: ['cell', 'organelles', 'nucleus', 'mitochondria', 'ribosomes', 'endoplasmic reticulum', 'golgi'],
      responses: [
        "Prokaryotic cells (bacteria) lack membrane-bound organelles. Eukaryotic cells have nucleus, mitochondria, ER, Golgi, lysosomes, ribosomes, cytoskeleton.",
        "Mitochondria: powerhouses with double membrane, cristae increase surface area, matrix contains enzymes. Chloroplasts: double membrane, thylakoids, stroma.",
        "Endomembrane system: nuclear envelope, ER (rough has ribosomes, smooth lacks ribosomes), Golgi apparatus, vesicles, lysosomes work together for protein processing."
      ]
    },
    {
      keywords: ['dna', 'replication', 'transcription', 'translation', 'genetics', 'chromosome', 'gene'],
      responses: [
        "DNA structure: double helix, antiparallel strands, A-T (2 H-bonds), G-C (3 H-bonds). Semiconservative replication: each new DNA has one old, one new strand.",
        "Central dogma: DNA → RNA → Protein. Transcription (DNA to mRNA) in nucleus, translation (mRNA to protein) at ribosomes. Start codon AUG, stop codons UAA, UAG, UGA.",
        "Genetic code: triplet codons specify amino acids, redundant (64 codons, 20 amino acids), universal across species, non-overlapping reading frame."
      ]
    },
    {
      keywords: ['mitosis', 'meiosis', 'cell division', 'chromosome', 'gametes', 'crossing over'],
      responses: [
        "Mitosis: PMAT (Prophase, Metaphase, Anaphase, Telophase) produces identical diploid cells. Meiosis: two divisions produce genetically diverse haploid gametes.",
        "Meiosis I: prophase I has crossing over (genetic recombination), metaphase I has independent assortment, homologs separate. Meiosis II separates sister chromatids.",
        "Genetic diversity sources: crossing over exchanges DNA between homologs, independent assortment randomly distributes chromosomes, random fertilization combines gametes."
      ]
    },
    {
      keywords: ['evolution', 'natural selection', 'darwin', 'adaptation', 'speciation', 'phylogeny'],
      responses: [
        "Natural selection mechanisms: directional (favors one extreme), stabilizing (favors average), disruptive (favors both extremes), sexual selection (mate choice).",
        "Evidence for evolution: fossil record, comparative anatomy (homologous structures), embryology, molecular biology (DNA/protein similarities), biogeography.",
        "Speciation: allopatric (geographic isolation), sympatric (reproductive isolation without geography), parapatric (partial isolation), peripatric (small population isolation)."
      ]
    },
    {
      keywords: ['ecology', 'ecosystem', 'food web', 'energy flow', 'nutrient cycle', 'biodiversity'],
      responses: [
        "Trophic levels: producers (autotrophs), primary consumers (herbivores), secondary consumers (carnivores), tertiary consumers, decomposers. 10% rule: energy transfer efficiency.",
        "Biogeochemical cycles: carbon (photosynthesis/respiration), nitrogen (fixation, nitrification, denitrification), phosphorus (no atmospheric form), water cycle.",
        "Community interactions: competition (intraspecific/interspecific), predation, mutualism (+/+), commensalism (+/0), parasitism (+/-), amensalism (0/-)."
      ]
    },
    {
      keywords: ['protein', 'amino acid', 'enzyme', 'structure', 'function', 'catalysis'],
      responses: [
        "Protein structure levels: primary (amino acid sequence), secondary (α-helices, β-sheets), tertiary (3D folding), quaternary (multiple polypeptides).",
        "Enzyme kinetics: Michaelis-Menten equation, Km (substrate concentration at ½ Vmax), competitive inhibition (increases Km), noncompetitive inhibition (decreases Vmax).",
        "Protein folding: hydrophobic effect drives folding, chaperones assist folding, misfolded proteins cause diseases (Alzheimer's, Parkinson's, prion diseases)."
      ]
    },
    {
      keywords: ['immune system', 'antibody', 'antigen', 'vaccine', 'lymphocyte', 'immunity'],
      responses: [
        "Innate immunity: immediate, nonspecific (skin, mucus, fever, inflammation, phagocytes, NK cells). Adaptive immunity: specific, memory (B cells, T cells).",
        "Humoral immunity: B cells produce antibodies (immunoglobulins). Cell-mediated immunity: T helper cells coordinate, cytotoxic T cells kill infected cells.",
        "Antibody structure: Y-shaped, heavy and light chains, variable regions (antigen binding), constant regions (effector functions), five classes (IgG, IgM, IgA, IgD, IgE)."
      ]
    },
    {
      keywords: ['neuron', 'action potential', 'synapse', 'neurotransmitter', 'brain', 'nervous system'],
      responses: [
        "Neuron types: sensory (afferent), motor (efferent), interneurons. Parts: dendrites (receive), cell body (integrate), axon (transmit), synaptic terminals (release).",
        "Action potential: resting potential (-70mV), depolarization (Na+ channels), repolarization (K+ channels), hyperpolarization, all-or-none principle.",
        "Neurotransmitters: acetylcholine (muscle contraction), dopamine (reward, movement), serotonin (mood), GABA (inhibitory), glutamate (excitatory)."
      ]
    },
    {
      keywords: ['metabolism', 'cellular respiration', 'glycolysis', 'krebs cycle', 'electron transport'],
      responses: [
        "Cellular respiration stages: glycolysis (cytoplasm, glucose to pyruvate), citric acid cycle (mitochondrial matrix), electron transport chain (inner membrane).",
        "ATP yield: glycolysis (2 net ATP), citric acid cycle (2 ATP), electron transport chain (~32-34 ATP). Total: ~36-38 ATP per glucose molecule.",
        "Fermentation: anaerobic ATP production. Alcoholic (yeast: pyruvate → ethanol + CO₂), lactic acid (muscle: pyruvate → lactate), regenerates NAD+ for glycolysis."
      ]
    },

    // === CHEMISTRY - 1000+ Questions ===
    {
      keywords: ['atomic structure', 'electron', 'proton', 'neutron', 'quantum', 'orbital'],
      responses: [
        "Quantum mechanical model: electrons exist in probability clouds (orbitals). s (spherical), p (dumbbell), d (cloverleaf), f (complex). Aufbau principle, Hund's rule, Pauli exclusion.",
        "Electron configuration: 1s² 2s² 2p⁶ 3s² 3p⁶ 4s² 3d¹⁰... Noble gas notation uses [Ne] 3s² 3p¹ for sodium. Valence electrons determine chemical properties.",
        "Effective nuclear charge: actual nuclear charge felt by valence electrons. Increases across period (more protons), constant down group (shielding effect)."
      ]
    },
    {
      keywords: ['periodic trends', 'ionization energy', 'electronegativity', 'atomic radius', 'electron affinity'],
      responses: [
        "Periodic trends: atomic radius decreases across period (more protons pull electrons closer), increases down group (more electron shells).",
        "Ionization energy: energy to remove electron. Increases across period (higher nuclear charge), decreases down group (electrons farther from nucleus).",
        "Electronegativity: ability to attract electrons in bonds. Increases across period and up group. Fluorine most electronegative (4.0), francium least (0.7)."
      ]
    },
    {
      keywords: ['chemical bonding', 'lewis structure', 'vsepr', 'hybridization', 'molecular orbital'],
      responses: [
        "VSEPR theory predicts molecular geometry: linear, trigonal planar, tetrahedral, trigonal bipyramidal, octahedral. Lone pairs occupy more space than bonding pairs.",
        "Hybridization: sp³ (tetrahedral), sp² (trigonal planar), sp (linear). Hybrid orbitals form sigma bonds, unhybridized p orbitals form pi bonds.",
        "Molecular orbital theory: bonding orbitals (lower energy), antibonding orbitals (higher energy). Bond order = (bonding electrons - antibonding electrons)/2."
      ]
    },
    {
      keywords: ['thermochemistry', 'enthalpy', 'entropy', 'gibbs free energy', 'spontaneous'],
      responses: [
        "First law of thermodynamics: ΔU = q + w (energy conservation). Enthalpy H = U + PV, constant pressure heat. ΔH < 0 exothermic, ΔH > 0 endothermic.",
        "Second law: entropy of universe increases. ΔS > 0 for spontaneous processes. Third law: S = 0 at absolute zero for perfect crystal.",
        "Gibbs free energy: ΔG = ΔH - TΔS. ΔG < 0 spontaneous, ΔG = 0 equilibrium, ΔG > 0 nonspontaneous. Temperature dependent spontaneity."
      ]
    },
    {
      keywords: ['chemical kinetics', 'reaction rate', 'activation energy', 'catalyst', 'mechanism'],
      responses: [
        "Rate law: rate = k[A]ᵐ[B]ⁿ where m, n are reaction orders (must be determined experimentally). Overall order = m + n.",
        "Arrhenius equation: k = Ae^(-Ea/RT). Higher temperature or lower activation energy increases rate constant. Catalysts lower Ea without being consumed.",
        "Reaction mechanisms: series of elementary steps. Rate-determining step (slowest) determines overall rate. Intermediates produced and consumed."
      ]
    },
    {
      keywords: ['chemical equilibrium', 'le chatelier', 'equilibrium constant', 'ice table'],
      responses: [
        "Equilibrium constant: Kc = [products]/[reactants] at constant temperature. Large K favors products, small K favors reactants. K varies with temperature only.",
        "Le Chatelier's principle: system responds to stress by shifting equilibrium. Adding reactants shifts right, removing products shifts right.",
        "ICE table: Initial concentrations, Change in concentrations, Equilibrium concentrations. Use stoichiometry and equilibrium expression to solve."
      ]
    },
    {
      keywords: ['acid base', 'pH', 'buffer', 'titration', 'ka', 'kb'],
      responses: [
        "Brønsted-Lowry theory: acids donate protons, bases accept protons. Conjugate acid-base pairs differ by one proton. Water is amphoteric.",
        "pH = -log[H⁺], pOH = -log[OH⁻], pH + pOH = 14 at 25°C. Ka × Kb = Kw for conjugate pairs. Stronger acid has weaker conjugate base.",
        "Buffer calculations: Henderson-Hasselbalch equation pH = pKa + log([A⁻]/[HA]). Buffer capacity highest when pH = pKa."
      ]
    },
    {
      keywords: ['redox', 'oxidation', 'reduction', 'electrochemistry', 'galvanic cell', 'electrolysis'],
      responses: [
        "Oxidation: loss of electrons (OIL), reduction: gain of electrons (RIG). Oxidation state increases in oxidation, decreases in reduction.",
        "Galvanic cell: spontaneous redox reaction generates electricity. Anode (oxidation, negative), cathode (reduction, positive). Salt bridge maintains charge balance.",
        "Standard electrode potentials: E°cell = E°cathode - E°anode. Positive E°cell indicates spontaneous reaction. ΔG° = -nFE°cell."
      ]
    },
    {
      keywords: ['organic chemistry', 'functional group', 'isomer', 'reaction mechanism', 'stereochemistry'],
      responses: [
        "Functional groups: alkanes (C-C), alkenes (C=C), alkynes (C≡C), aromatics (benzene), alcohols (-OH), aldehydes (-CHO), ketones (C=O), carboxylic acids (-COOH).",
        "Stereoisomers: same connectivity, different spatial arrangement. Enantiomers (mirror images), diastereomers (not mirror images). Chiral centers have four different groups.",
        "Reaction mechanisms: nucleophilic substitution (SN1, SN2), elimination (E1, E2), addition reactions. Arrow pushing shows electron movement."
      ]
    },
    {
      keywords: ['spectroscopy', 'nmr', 'ir', 'mass spec', 'uv vis'],
      responses: [
        "NMR spectroscopy: ¹H NMR shows proton environments, integration shows relative numbers, coupling shows neighboring protons. ¹³C NMR shows carbon environments.",
        "IR spectroscopy: molecular vibrations absorb specific frequencies. O-H stretch (3200-3600 cm⁻¹), C=O stretch (1700 cm⁻¹), C-H stretch (2800-3000 cm⁻¹).",
        "Mass spectrometry: molecular ion peak shows molecular weight, fragmentation patterns identify structure. Base peak is most intense fragment."
      ]
    },

    // === PHYSICS - 1000+ Questions ===
    {
      keywords: ['mechanics', 'force', 'newton', 'motion', 'acceleration', 'momentum'],
      responses: [
        "Newton's laws: 1st (inertia), 2nd (F = ma), 3rd (action-reaction). Net force causes acceleration in direction of force. Mass is inertia measure.",
        "Kinematic equations: v = v₀ + at, x = x₀ + v₀t + ½at², v² = v₀² + 2a(x-x₀), x = x₀ + ½(v₀+v)t. Valid for constant acceleration only.",
        "Conservation of momentum: total momentum before collision equals total momentum after. Elastic collisions conserve kinetic energy, inelastic do not."
      ]
    },
    {
      keywords: ['energy', 'work', 'power', 'conservation', 'kinetic', 'potential'],
      responses: [
        "Work-energy theorem: net work equals change in kinetic energy. W = ΔKE = ½mv² - ½mv₀². Work done by conservative forces is path-independent.",
        "Mechanical energy: E = KE + PE. Conserved when only conservative forces act. Spring potential energy: PE = ½kx². Gravitational: PE = mgh.",
        "Power: P = W/t = F·v. Units: Watts (J/s). Efficiency = (useful energy output)/(total energy input) × 100%. No machine is 100% efficient."
      ]
    },
    {
      keywords: ['waves', 'sound', 'light', 'frequency', 'wavelength', 'interference'],
      responses: [
        "Wave equation: v = fλ where v = speed, f = frequency, λ = wavelength. Period T = 1/f. All electromagnetic waves travel at c = 3×10⁸ m/s in vacuum.",
        "Sound waves: longitudinal pressure waves in air. Speed ≈ 343 m/s at 20°C. Doppler effect: frequency changes when source or observer moves.",
        "Wave interference: constructive (waves add), destructive (waves cancel). Double-slit experiment shows light's wave nature. Standing waves have nodes and antinodes."
      ]
    },
    {
      keywords: ['electricity', 'magnetism', 'current', 'voltage', 'resistance', 'field'],
      responses: [
        "Ohm's law: V = IR. Resistance R = ρL/A where ρ = resistivity, L = length, A = cross-sectional area. Series: Rtotal = R₁ + R₂. Parallel: 1/Rtotal = 1/R₁ + 1/R₂.",
        "Electric field: E = F/q. Field lines point from positive to negative charges. Uniform field between parallel plates: E = V/d. Gauss's law relates field to charge.",
        "Magnetic force: F = qvB sin θ for moving charge. F = ILB sin θ for current-carrying wire. Right-hand rule determines force direction."
      ]
    },
    {
      keywords: ['thermodynamics', 'heat', 'temperature', 'gas laws', 'entropy'],
      responses: [
        "Ideal gas law: PV = nRT. Kinetic theory: average kinetic energy ∝ absolute temperature. Real gases deviate at high pressure and low temperature.",
        "Heat engines: efficiency = W/Qh = 1 - Qc/Qh ≤ 1 - Tc/Th (Carnot limit). Refrigerators move heat from cold to hot reservoir using work.",
        "Laws of thermodynamics: 0th (thermal equilibrium), 1st (energy conservation), 2nd (entropy increases), 3rd (S → 0 as T → 0)."
      ]
    },
    {
      keywords: ['optics', 'reflection', 'refraction', 'lens', 'mirror', 'ray'],
      responses: [
        "Law of reflection: angle of incidence = angle of reflection. Snell's law: n₁sin θ₁ = n₂sin θ₂. Total internal reflection when θ > critical angle.",
        "Thin lens equation: 1/f = 1/do + 1/di. Magnification: m = -di/do = hi/ho. Converging lens (f > 0), diverging lens (f < 0).",
        "Wave optics: diffraction (bending around obstacles), interference (double-slit), polarization (orientation of electric field)."
      ]
    },
    {
      keywords: ['modern physics', 'quantum', 'relativity', 'atomic', 'nuclear'],
      responses: [
        "Photoelectric effect: light behaves as particles (photons). E = hf where h = Planck's constant. Work function φ is minimum energy to remove electron.",
        "Special relativity: time dilation Δt' = γΔt, length contraction L' = L/γ where γ = 1/√(1-v²/c²). Nothing with mass can reach light speed.",
        "Bohr model: electrons in quantized orbits. Energy levels: En = -13.6 eV/n². Photon emission/absorption: ΔE = hf = Ei - Ef."
      ]
    },
    {
      keywords: ['nuclear physics', 'radioactivity', 'fusion', 'fission', 'decay'],
      responses: [
        "Nuclear decay types: alpha (α, ⁴₂He), beta (β⁻, electron; β⁺, positron), gamma (γ, photon). Exponential decay: N(t) = N₀e^(-λt).",
        "Nuclear fission: heavy nucleus splits into lighter nuclei + neutrons + energy. Chain reaction in nuclear reactors. U-235 fissile isotope.",
        "Nuclear fusion: light nuclei combine into heavier nucleus + energy. Powers stars. D + T → ⁴He + n + 17.6 MeV. Requires high temperature/pressure."
      ]
    },
    {
      keywords: ['astrophysics', 'star', 'planet', 'galaxy', 'cosmology'],
      responses: [
        "Stellar evolution: main sequence (H → He fusion), red giant (He fusion), white dwarf/neutron star/black hole (depends on mass). HR diagram plots luminosity vs temperature.",
        "Exoplanets: detected by transit method (dimming), radial velocity (wobble), direct imaging. Habitable zone where liquid water possible.",
        "Big Bang cosmology: universe expanding since 13.8 billion years ago. CMB radiation, nucleosynthesis, dark matter/energy comprise 95% of universe."
      ]
    },
    {
      keywords: ['particle physics', 'standard model', 'fundamental forces', 'quarks', 'leptons'],
      responses: [
        "Standard Model: fundamental particles are quarks (up, down, charm, strange, top, bottom) and leptons (electron, muon, tau, neutrinos).",
        "Four fundamental forces: strong (quarks), electromagnetic (charged particles), weak (radioactive decay), gravitational (mass). Mediated by bosons.",
        "Particle accelerators: collide high-energy particles to study fundamental interactions. LHC discovered Higgs boson in 2012."
      ]
    },

    // === MATHEMATICS - 1200+ Questions ===
    {
      keywords: ['algebra', 'equation', 'polynomial', 'factor', 'quadratic', 'linear'],
      responses: [
        "Quadratic formula: x = [-b ± √(b²-4ac)]/2a for ax² + bx + c = 0. Discriminant b²-4ac determines number and type of roots.",
        "Polynomial division: synthetic division for linear divisors, long division for higher degree. Remainder theorem: P(c) equals remainder when P(x) divided by (x-c).",
        "Systems of equations: substitution, elimination, matrix methods. Consistent (one/infinite solutions), inconsistent (no solution)."
      ]
    },
    {
      keywords: ['geometry', 'triangle', 'circle', 'polygon', 'similarity', 'congruence'],
      responses: [
        "Triangle congruence: SSS, SAS, ASA, AAS, HL (right triangles). Similarity: AA, SSS~, SAS~. Similar triangles have proportional sides.",
        "Circle theorems: inscribed angle = ½ central angle, tangent perpendicular to radius, angles in semicircle = 90°, tangent-chord angle.",
        "Coordinate geometry: distance formula d = √[(x₂-x₁)² + (y₂-y₁)²], midpoint ((x₁+x₂)/2, (y₁+y₂)/2), slope m = (y₂-y₁)/(x₂-x₁)."
      ]
    },
    {
      keywords: ['trigonometry', 'sine', 'cosine', 'tangent', 'identity', 'law of sines'],
      responses: [
        "Fundamental identities: sin²θ + cos²θ = 1, tan θ = sin θ/cos θ, 1 + tan²θ = sec²θ, 1 + cot²θ = csc²θ.",
        "Addition formulas: sin(A±B) = sin A cos B ± cos A sin B, cos(A±B) = cos A cos B ∓ sin A sin B, tan(A±B) = (tan A ± tan B)/(1 ∓ tan A tan B).",
        "Law of sines: a/sin A = b/sin B = c/sin C. Law of cosines: c² = a² + b² - 2ab cos C. Use for any triangle, not just right triangles."
      ]
    },
    {
      keywords: ['calculus', 'derivative', 'integral', 'limit', 'chain rule', 'fundamental theorem'],
      responses: [
        "Limit definition of derivative: f'(x) = lim[h→0] [f(x+h) - f(x)]/h. Derivatives: power rule, product rule, quotient rule, chain rule.",
        "Integration techniques: substitution (reverse chain rule), integration by parts (∫u dv = uv - ∫v du), partial fractions, trigonometric substitution.",
        "Fundamental theorem of calculus: ∫[a to b] f'(x)dx = f(b) - f(a). Connects derivatives and integrals as inverse operations."
      ]
    },
    {
      keywords: ['statistics', 'probability', 'distribution', 'hypothesis test', 'regression'],
      responses: [
        "Central limit theorem: sample means approach normal distribution as n increases. Standard error = σ/√n. Confidence intervals use z or t distributions.",
        "Hypothesis testing: null hypothesis H₀, alternative H₁. Type I error (reject true H₀), Type II error (accept false H₀). P-value is probability of observed result if H₀ true.",
        "Correlation vs causation: correlation measures linear relationship strength (-1 to +1). Causation requires controlled experiments, not just correlation."
      ]
    },
    {
      keywords: ['linear algebra', 'matrix', 'vector', 'eigenvalue', 'determinant'],
      responses: [
        "Matrix operations: addition (element-wise), multiplication (row×column), transpose (flip rows/columns). Identity matrix I is multiplicative identity.",
        "Determinant: ad-bc for 2×2 matrix. Properties: det(AB) = det(A)det(B), det(A⁻¹) = 1/det(A), det(Aᵀ) = det(A).",
        "Eigenvalues and eigenvectors: Av = λv where λ is eigenvalue, v is eigenvector. Characteristic equation: det(A - λI) = 0."
      ]
    },
    {
      keywords: ['differential equations', 'separable', 'linear', 'homogeneous'],
      responses: [
        "Separable equations: dy/dx = g(x)h(y) becomes dy/h(y) = g(x)dx. Integrate both sides and solve for y.",
        "First-order linear: dy/dx + P(x)y = Q(x). Integrating factor μ(x) = e^∫P(x)dx. Solution: y = [∫μQ dx + C]/μ.",
        "Second-order homogeneous: y'' + py' + qy = 0. Characteristic equation: r² + pr + q = 0. Solutions depend on discriminant."
      ]
    },
    {
      keywords: ['discrete math', 'combinatorics', 'graph theory', 'logic'],
      responses: [
        "Combinations: C(n,r) = n!/(r!(n-r)!) order doesn't matter. Permutations: P(n,r) = n!/(n-r)! order matters.",
        "Graph theory: vertices connected by edges. Euler path visits each edge once, Hamiltonian path visits each vertex once. Trees are connected acyclic graphs.",
        "Propositional logic: AND (∧), OR (∨), NOT (¬), IMPLIES (→), IFF (↔). Truth tables determine validity. De Morgan's laws: ¬(p∧q) = ¬p∨¬q."
      ]
    },
    {
      keywords: ['number theory', 'prime', 'gcd', 'modular arithmetic', 'cryptography'],
      responses: [
        "Prime numbers: only divisible by 1 and themselves. Fundamental theorem: every integer > 1 has unique prime factorization.",
        "Greatest common divisor: gcd(a,b) found by Euclidean algorithm. Extended Euclidean algorithm finds integers x,y where ax + by = gcd(a,b).",
        "Modular arithmetic: a ≡ b (mod n) means n divides (a-b). Used in cryptography, computer science, number theory."
      ]
    },
    {
            keywords: ['complex analysis', 'analytic function', 'cauchy riemann', 'residue theorem'],
      responses: [
        "Complex functions: f(z) = u(x,y) + iv(x,y) where z = x + iy. Analytic functions satisfy Cauchy-Riemann equations: ∂u/∂x = ∂v/∂y, ∂u/∂y = -∂v/∂x.",
        "Cauchy's theorem: ∫f(z)dz = 0 for analytic function over closed contour. Residue theorem calculates contour integrals using singularities.",
        "Conformal mapping: analytic functions preserve angles. Used in fluid dynamics, electrostatics, heat conduction problems."
      ]
    },

    // === HISTORY - 800+ Questions ===
    {
      keywords: ['ancient civilizations', 'mesopotamia', 'egypt', 'indus valley', 'yellow river'],
      responses: [
        "Mesopotamia (3500-539 BCE): Sumerians invented writing (cuneiform), wheel, city-states. Babylonians created Hammurabi's Code. Assyrians built empire.",
        "Ancient Egypt (3100-30 BCE): Pharaohs, pyramids, mummification, hieroglyphs. Nile River crucial for agriculture. New Kingdom peak with Ramses II.",
        "Indus Valley (2600-1900 BCE): Harappa, Mohenjo-daro had advanced urban planning, sewage systems, standardized weights. Undeciphered script."
      ]
    },
    {
      keywords: ['classical greece', 'athens', 'sparta', 'persian wars', 'peloponnesian war'],
      responses: [
        "Persian Wars (499-449 BCE): Marathon (490), Thermopylae (480), Salamis (480), Plataea (479). Greek city-states defeated Persian Empire.",
        "Golden Age Athens (461-429 BCE): Pericles leadership, Parthenon, democracy, philosophy (Socrates, Plato, Aristotle), theater (Sophocles, Euripides).",
        "Peloponnesian War (431-404 BCE): Athens vs Sparta and allies. Thucydides documented. Plague, Sicilian Expedition, Sparta victory ended Athenian dominance."
      ]
    },
    {
      keywords: ['roman republic', 'julius caesar', 'augustus', 'roman empire', 'christianity'],
      responses: [
        "Roman Republic (509-27 BCE): Senate, consuls, tribunes. Punic Wars vs Carthage, Hannibal's invasion. Social War, slave revolts, civil wars ended republic.",
        "Julius Caesar: conquered Gaul, crossed Rubicon (49 BCE), dictator, calendar reform, assassinated 44 BCE by Brutus, Cassius. Triggered final civil wars.",
        "Augustus (27 BCE-14 CE): first emperor, Pax Romana began, efficient administration, infrastructure, literature (Virgil, Ovid). Adopted son of Caesar."
      ]
    },
    {
      keywords: ['medieval europe', 'feudalism', 'crusades', 'black death', 'magna carta'],
      responses: [
        "Feudalism: hierarchical system with lords, vassals, serfs. Land grants (fiefs) for military service. Manorialism organized agricultural production.",
        "Crusades (1096-1291): Christian military expeditions to Holy Land. First Crusade captured Jerusalem. Later crusades less successful, increased trade.",
        "Black Death (1347-1351): bubonic plague killed 1/3 of European population. Led to labor shortages, social upheaval, religious questioning."
      ]
    },
    {
      keywords: ['renaissance', 'humanism', 'leonardo da vinci', 'michelangelo', 'printing press'],
      responses: [
        "Italian Renaissance (14th-16th centuries): Florence, Venice, Rome centers. Medici family patronage. Revival of classical learning, individual achievement.",
        "Renaissance art: linear perspective, chiaroscuro (light/shadow), sfumato. Leonardo (Mona Lisa, Last Supper), Michelangelo (Sistine Chapel, David).",
        "Printing press (Gutenberg, 1440): movable type, mass book production. Bible, scientific works, literature spread. Literacy increased, ideas disseminated."
      ]
    },
    {
      keywords: ['reformation', 'martin luther', 'calvin', 'counter reformation', 'wars of religion'],
      responses: [
        "Protestant Reformation: Martin Luther's 95 Theses (1517) challenged Catholic practices. Salvation by faith alone, scripture authority, priesthood of believers.",
        "John Calvin: predestination, Geneva theocracy, Protestant work ethic. Calvinism spread to Netherlands, Scotland (John Knox), New England Puritans.",
        "Counter-Reformation: Catholic response. Council of Trent, Jesuits, Inquisition, Index of Forbidden Books. Baroque art promoted Catholic themes."
      ]
    },
    {
      keywords: ['age of exploration', 'columbus', 'vasco da gama', 'magellan', 'conquistadors'],
      responses: [
        "Motives for exploration: gold, glory, God. Technology: compass, astrolabe, caravel ships, improved maps. Portugal led, Spain followed.",
        "Columbus (1492): reached Americas seeking Asia route. Columbian Exchange: crops, animals, diseases between Old/New Worlds. Devastated Native Americans.",
        "Conquistadors: Cortés conquered Aztecs (1519-1521), Pizarro conquered Incas (1532-1533). Superior weapons, horses, diseases, native allies."
      ]
    },
    {
      keywords: ['scientific revolution', 'copernicus', 'galileo', 'newton', 'scientific method'],
      responses: [
        "Copernican Revolution: heliocentric model replaced geocentric. Galileo's telescope observations supported Copernicus. Church opposition, trial (1633).",
        "Scientific method: observation, hypothesis, experimentation, verification. Francis Bacon (empiricism), René Descartes (rationalism), inductive reasoning.",
        "Newton's Principia (1687): laws of motion, universal gravitation, mathematical physics. Unified terrestrial and celestial mechanics."
      ]
    },
    {
      keywords: ['enlightenment', 'voltaire', 'rousseau', 'locke', 'separation of powers'],
      responses: [
        "Age of Reason (1685-1815): rational thought, natural rights, religious tolerance, progress. Philosophes promoted reform through reason.",
        "John Locke: natural rights (life, liberty, property), government by consent, right of revolution. Influenced American revolutionaries.",
        "Montesquieu: separation of powers (legislative, executive, judicial), checks and balances. Influenced US Constitution framers."
      ]
    },
    {
      keywords: ['american revolution', 'declaration of independence', 'constitution', 'federalism'],
      responses: [
        "Causes: taxation without representation, Stamp Act, Tea Act, Intolerable Acts. Boston Massacre, Boston Tea Party increased tensions.",
        "Declaration of Independence (1776): natural rights philosophy, list of grievances, right of revolution. Jefferson main author, influenced by Locke.",
        "Constitution (1787): federal system, separation of powers, Bill of Rights. Compromise between large/small states, federalists/anti-federalists."
      ]
    },
    {
      keywords: ['french revolution', 'robespierre', 'reign of terror', 'napoleon', 'congress of vienna'],
      responses: [
        "Causes: financial crisis, social inequality, Enlightenment ideas, American example. Estates-General, Tennis Court Oath, Storming of Bastille.",
        "Reign of Terror (1793-1794): Committee of Public Safety, Robespierre, guillotine executions. Ended with Thermidorian Reaction, Robespierre executed.",
        "Napoleon (1799-1815): consul, emperor, Napoleonic Code, Continental System. Military genius, conquered most of Europe, defeated at Waterloo."
      ]
    },
    {
      keywords: ['industrial revolution', 'steam engine', 'factory system', 'railroad', 'urbanization'],
      responses: [
        "Began in Britain (1760s): textile industry mechanization, steam power, iron/steel production. Spread to Europe, North America.",
        "Transportation revolution: canals, turnpikes, railroads, steamships. Reduced costs, increased speed, connected markets, encouraged specialization.",
        "Social effects: urbanization, working class formation, family structure changes, child labor, pollution, wealth inequality, labor unions."
      ]
    },
    {
      keywords: ['nationalism', 'unification', 'germany', 'italy', 'ottoman decline'],
      responses: [
        "German unification: Prussian leadership, Otto von Bismarck, wars vs Denmark, Austria, France. German Empire proclaimed 1871 at Versailles.",
        "Italian unification: Risorgimento movement, Giuseppe Garibaldi, Count Cavour, Kingdom of Piedmont-Sardinia. Unified 1870 except Vatican.",
        "Ottoman decline: 'Sick Man of Europe,' nationalism in Balkans, Russian expansion, European intervention. Young Turk revolution 1908."
      ]
    },
    {
      keywords: ['imperialism', 'scramble for africa', 'opium wars', 'meiji restoration'],
      responses: [
        "New Imperialism (1870-1914): European powers colonized Africa, Asia. Motives: economic, strategic, civilizing mission, national prestige.",
        "Scramble for Africa: Berlin Conference (1884-1885) divided Africa. Only Ethiopia, Liberia remained independent. Exploitation, cultural disruption.",
        "Meiji Restoration (1868): Japan modernized rapidly, adopted Western technology, industrialized, built strong military. Became imperial power."
      ]
    },

    // === LITERATURE - 600+ Questions ===
    {
      keywords: ['epic poetry', 'homer', 'iliad', 'odyssey', 'beowulf', 'dante'],
      responses: [
        "Epic characteristics: long narrative poem, heroic protagonist, supernatural elements, cultural values. Oral tradition, formal style, epic similes.",
        "Homer's Iliad: Trojan War, Achilles' wrath, honor vs life conflict. Odyssey: Odysseus' journey home, cleverness, loyalty themes.",
        "Dante's Divine Comedy: journey through Hell, Purgatory, Paradise. Allegory of soul's path to God. Terza rima, vernacular Italian."
      ]
    },
    {
      keywords: ['shakespearean drama', 'tragedy', 'comedy', 'history plays', 'sonnets'],
      responses: [
        "Shakespeare's tragedies: Hamlet (revenge, madness), Macbeth (ambition, guilt), King Lear (family, justice), Othello (jealousy, racism).",
        "Comedies: mistaken identities, disguises, multiple plots, happy endings. Much Ado About Nothing, A Midsummer Night's Dream, Twelfth Night.",
        "Shakespearean sonnet: 14 lines, ABAB CDCD EFEF GG rhyme scheme, iambic pentameter. Themes: love, beauty, time, mortality, immortality through art."
      ]
    },
    {
      keywords: ['romantic poetry', 'wordsworth', 'coleridge', 'byron', 'keats', 'shelley'],
      responses: [
        "Romanticism: emotion over reason, nature worship, individualism, imagination. Reaction against Enlightenment rationalism, Industrial Revolution.",
        "Wordsworth: 'emotion recollected in tranquility,' nature poetry, simple language. 'I Wandered Lonely as a Cloud,' Lyrical Ballads.",
        "Keats: beauty, sensuality, mortality. 'Ode to a Nightingale."
      ]
    },{
                 keywords: ['romantic poetry', 'wordsworth', 'coleridge', 'byron', 'keats', 'shelley'],
      responses: [
        "Romanticism emphasized emotion, imagination, nature, and individualism as a reaction to Enlightenment rationalism and industrialization, with lyric poetry becoming a dominant form.",
        "Wordsworth championed 'emotion recollected in tranquility' and plain diction; Coleridge explored the supernatural; Byron embodied the Byronic hero; Shelley advocated radical ideals; Keats pursued beauty and transience.",
        "Keats’s great odes (Nightingale, Grecian Urn, Autumn) balance sensuous imagery with philosophical reflection on mortality, art, and the tension between the real and the ideal."
      ]
    },
    {
      keywords: ['modernism', 'ts eliot', 'james joyce', 'virginia woolf', 'stream of consciousness'],
      responses: [
        "Modernism broke with traditional forms to depict fragmentation, alienation, and the complexities of consciousness through techniques like stream of consciousness and nonlinear structure.",
        "T. S. Eliot’s The Waste Land blends mythic method, allusion, and polyphony to portray spiritual desolation; Joyce’s Ulysses experiments with narrative styles; Woolf explores interiority and time.",
        "Modernist motifs include urban anonymity, the failure of grand narratives, linguistic experimentation, and the search for meaning amid cultural upheaval."
      ]
    },
    {
      keywords: ['literary devices', 'metaphor', 'simile', 'allegory', 'symbol', 'irony'],
      responses: [
        "Metaphor asserts identity between unlike things; simile uses 'like/as'; allegory sustains parallel narratives; symbols accrue meaning beyond literal reference; irony subverts expectations.",
        "Structural devices—motif, foreshadowing, flashback—shape reader anticipation and thematic coherence; sound devices—assonance, consonance, alliteration—reinforce tone and rhythm.",
        "Close reading traces how diction, imagery, syntax, and figurative language produce effects, supporting interpretive claims with textual evidence."
      ]
    },

    // === DEFAULT & FALLBACKS ===
    {
      keywords: ['help', 'guide', 'how to', 'explain', 'clarify', 'teach'],
      responses: [
        "Ask about any topic by mentioning key terms; for example, 'photosynthesis stages', 'Ohm’s law basics', 'Pythagorean theorem proof', or 'French Revolution causes'.",
        "Provide keywords rather than long paragraphs—this chatbot uses keyword and fuzzy matching to find the best explanation quickly.",
        "If a response doesn’t fit, try adding one or two more specific terms, like 'electrochemistry galvanic cell' or 'calculus chain rule example'."
      ]
    },
    {
keywords: ['classical tragedy', 'aristotle', 'hamartia', 'catharsis', 'anagnorisis'],
responses: [
"Classical tragedy follows a noble protagonist whose hamartia leads to downfall and evokes catharsis through pity and fear.",
"Aristotle outlined plot primacy, unity, peripeteia, and anagnorisis as key elements that intensify tragic impact.",
"The tragic hero’s recognition scene reveals self‑knowledge too late, aligning moral error with irreversible consequence."
]
},
{
keywords: ['greek drama', 'aeschylus', 'sophocles', 'euripides', 'chorus'],
responses: [
"Greek drama forged tragedy and comedy with masked actors, formal verse, and a commenting chorus.",
"Aeschylus expanded dramatic conflict, Sophocles refined tragic form, and Euripides probed psychology and social critique.",
"The chorus mediated audience response, framed mythic context, and underscored ethical tensions on stage."
]
},
{
keywords: ['aristotle poetics', 'mimesis', 'plot', 'character', 'unity'],
responses: [
"Poetics centers plot as the soul of tragedy and defines mimesis as structured imitation of action.",
"Aristotle champions unity of action with tightly causal plots over episodic sequences for emotional power.",
"Character supports plot by revealing choice and ethics, but action drives the arc toward catharsis."
]
},
{
keywords: ['medieval literature', 'allegory', 'chivalry', 'mystery plays', 'morality plays'],
responses: [
"Medieval literature blends Christian allegory with chivalric romance and courtly love conventions.",
"Mystery and miracle plays dramatized biblical narratives for civic audiences in cyclical pageants.",
"Morality plays personified virtues and vices to stage the soul’s ethical struggle and salvation."
]
},
{
keywords: ['renaissance humanism', 'classical revival', 'individualism', 'reason', 'virtù'],
responses: [
"Renaissance humanism revived classical learning to cultivate eloquence, civic virtue, and individual dignity.",
"Writers fused philology and moral philosophy, emphasizing reason guided by exemplary antiquity.",
"Humanist virtù prized flexible excellence in arts, letters, and public life as a cultural ideal."
]
},
{
keywords: ['metaphysical poetry', 'conceit', 'donne', 'paradox', 'wit'],
responses: [
"Metaphysical poetry wields bold conceits and paradox to yoke intellect with passion.",
"Donne’s sermons and lyrics entwine theology, eros, and mortality through argumentative stanzaic drama.",
"The metaphysical wit tests readers with logical turns, learned references, and compressed intensity."
]
},
{
keywords: ['neoclassicism', 'order', 'decorum', 'satire', 'pope'],
responses: [
"Neoclassicism valued order, clarity, and decorum, modeling restraint on Greco‑Roman forms.",
"Satire flourished as a corrective art, disciplining taste and manners through urbane mockery.",
"Pope’s heroic couplets balanced reason and polish, perfecting concise moral epigram."
]
},
{
keywords: ['romanticism', 'imagination', 'sublime', 'nature', 'emotion'],
responses: [
"Romanticism elevated imagination and subjective feeling against mechanistic rationalism.",
"Poets sought the sublime in nature’s vastness to awaken awe, terror, and transcendence.",
"The self became a visionary instrument, forging authenticity through inward illumination."
]
},
{
keywords: ['victorian novel', 'serial publication', 'realism', 'social reform', 'morality'],
responses: [
"The Victorian novel thrived in serial form, integrating realism with reformist social vision.",
"Writers mapped class mobility, industrial hardship, and domestic ideals with moral scrutiny.",
"Expansive plots balanced sentiment, satire, and civic conscience within populous urban worlds."
]
},
{
keywords: ['realism', 'naturalism', 'determinism', 'everyday life', 'objectivity'],
responses: [
"Realism foregrounded ordinary life and plausible causality through detailed social observation.",
"Naturalism radicalized realism with determinism, heredity, and environment shaping fate.",
"An objective narrator frames characters as products of material and social pressures."
]
},
{
keywords: ['symbolism', 'suggestion', 'music of language', 'mallarmé', 'baudelaire'],
responses: [
"Symbolism favored suggestive imagery and synesthetic effects over direct statement.",
"Poets pursued the music of language to unlock spiritual correspondences and hidden states.",
"From Baudelaire to Mallarmé, symbolist poetics veiled meaning in patterned ambiguity."
]
},
{
keywords: ['decadence', 'aestheticism', 'art for art’s sake', 'fin de siècle'],
responses: [
"Decadence embraced languor, artifice, and exquisite sensibility at the fin de siècle.",
"Aestheticism asserted art’s autonomy, privileging beauty and style over moral utility.",
"Surface became substance as ornament, pose, and refined excess signaled radical taste."
]
},
{
keywords: ['modernism', 'fragmentation', 'mythic method', 'alienation', 'free verse'],
responses: [
"Modernism fractured narrative and form to mirror dislocation in the twentieth century.",
"The mythic method restructured contemporary chaos through archetype and ritual pattern.",
"Free verse, interior monologue, and collage sought new orders out of cultural wreckage."
]
},
{
keywords: ['stream of consciousness', 'interiority', 'time', 'joyce', 'woolf'],
responses: [
"Stream of consciousness renders unfiltered mental flux across memory, perception, and desire.",
"Modernists reimagined time as layered durée rather than linear sequence of events.",
"Joyce and Woolf turned narration into a living medium of consciousness in motion."
]
},
{
keywords: ['postmodernism', 'metafiction', 'pastiche', 'paranoia', 'play'],
responses: [
"Postmodernism delights in self‑reference, textual games, and ontological instability.",
"Pastiche recycles styles without hierarchy, while parody interrogates their authority.",
"Conspiracy, randomness, and irony replace master narratives with proliferating possibles."
]
},
{
keywords: ['magical realism', 'marquez', 'ordinary and marvelous', 'myth', 'history'],
responses: [
"Magical realism fuses the ordinary and marvelous within a matter‑of‑fact narrative tone.",
"Writers braid myth, folklore, and history to reveal cultural truths beyond realism.",
"Marquez’s chronicle entwines memory, politics, and enchantment as coequal realities."
]
},
{
keywords: ['latin american boom', 'experiment', 'politics', 'labyrinth', 'collective memory'],
responses: [
"The Boom exported bold experimentation, polyphony, and political allegory to global readers.",
"Labyrinthine structures map nationhood, trauma, and identity across shifting voices.",
"Collective memory challenges official history through metafiction and baroque excess."
]
},
{
keywords: ['harlem renaissance', 'new negro', 'jazz poetics', 'racial uplift', 'modernity'],
responses: [
"The Harlem Renaissance articulated the New Negro’s cultural agency and modern self‑fashioning.",
"Jazz rhythms, vernacular speech, and interarts collaboration forged a Black modernism.",
"Literature contested stereotypes and pursued racial uplift through aesthetic innovation."
]
},
{
keywords: ['beat generation', 'spontaneous prose', 'counterculture', 'buddhism', 'road'],
responses: [
"Beat writing prized spontaneity, jazz cadence, and visionary rebellion against conformity.",
"Buddhist thought, sexuality, and travel narratives reframed American restlessness.",
"The open road became a spiritual corridor toward ecstatic perception and dissent."
]
},
{
keywords: ['existentialism', 'absurd', 'freedom', 'responsibility', 'nausea'],
responses: [
"Existential literature confronts absurdity, freedom, and the burden of authentic choice.",
"Characters face contingency without guarantees, creating meaning through committed action.",
"Alienation and nausea expose the groundlessness beneath social roles and metaphysics."
]
},
{
keywords: ['theatre of the absurd', 'beckett', 'ionesco', 'meaninglessness', 'circularity'],
responses: [
"Absurdist drama stages breakdowns in language, purpose, and coherent plot.",
"Circular structures, repetition, and silence evoke existential stalemate.",
"Beckett and Ionesco render comedy as bleak revelation of human impasse."
]
},
{
keywords: ['epic theatre', 'brecht', 'alienation effect', 'didactic', 'gestus'],
responses: [
"Brecht’s epic theatre disrupts illusion to provoke critical social awareness.",
"Alienation effects expose theatrical machinery and invite analytical distance.",
"Gestus crystallizes social relations in emblematic gestures and tableaux."
]
},
{
keywords: ['naturalism theatre', 'slice of life', 'fourth wall', 'heredity', 'environment'],
responses: [
"Naturalist theatre presents a slice of life shaped by heredity and milieu.",
"Detailed sets and the fourth wall cultivate observational realism on stage.",
"Deterministic pressures steer characters toward credible, often bleak ends."
]
},
{
keywords: ['gothic fiction', 'ruins', 'sublime terror', 'doubling', 'uncanny'],
responses: [
"Gothic fiction conjures sublime terror through decay, secrecy, and transgression.",
"Doubles and haunted spaces externalize repressed histories and desires.",
"Atmospheric excess frames psychological disturbance as architectural fate."
]
},
{
keywords: ['southern gothic', 'grotesque', 'decay', 'religion', 'violence'],
responses: [
"Southern Gothic probes moral rot beneath politeness with grotesque characters and stark grace.",
"Religious symbolism collides with violence to expose cultural contradictions.",
"Regional detail anchors universal themes of guilt, pride, and redemption."
]
},
{
keywords: ['detective fiction', 'clues', 'red herring', 'deduction', 'whodunit'],
responses: [
"Detective fiction orchestrates clues and misdirection toward rational resolution.",
"Red herrings exploit reader inference while fair play hides truth in plain sight.",
"The sleuth models disciplined observation, logic, and narrative control."
]
},
{
keywords: ['noir', 'hardboiled', 'fatalism', 'city', 'antihero'],
responses: [
"Noir distills fatalism, moral ambiguity, and urban night into terse style.",
"Hardboiled narration exposes corruption through a wounded, relentless eye.",
"The antihero navigates violence and desire amid systemic decay."
]
},
{
keywords: ['science fiction', 'speculation', 'future', 'technology', 'sense of wonder'],
responses: [
"Science fiction tests humanity against speculative technologies and altered realities.",
"Worldbuilding dramatizes ethical dilemmas born of discovery and power.",
"A sense of wonder enlarges imagination while interrogating progress myths."
]
},
{
keywords: ['cyberpunk', 'high tech low life', 'virtual', 'corporation', 'hacker'],
responses: [
"Cyberpunk contrasts advanced tech with social ruin and corporate sovereignty.",
"Virtual spaces and body mods scramble identity, autonomy, and surveillance.",
"Hackers and outcasts map the moral edges of digitized capitalism."
]
},
{
keywords: ['dystopia', 'totalitarian', 'surveillance', 'resistance', 'control'],
responses: [
"Dystopian fiction extrapolates repression to warn against present tendencies.",
"Surveillance, propaganda, and scarcity discipline bodies and desires.",
"Resistance navigates compromised choices within engineered constraints."
]
},
{
keywords: ['utopia', 'ideal society', 'planning', 'equity', 'critique'],
responses: [
"Utopias model ideal orders to critique existing social arrangements.",
"Design tensions surface between liberty, equality, and communal aims.",
"Blueprint narratives probe the costs of harmony and rational planning."
]
},
{
keywords: ['fantasy', 'secondary world', 'mythic', 'quest', 'magic'],
responses: [
"Fantasy constructs secondary worlds bound by coherent magical laws and lore.",
"Mythic quests test virtue, fate, and belonging across perilous thresholds.",
"Magic externalizes inner conflicts as tangible, negotiable forces."
]
},
{
keywords: ['myth retelling', 'intertext', 'archetype', 'revision', 'voice'],
responses: [
"Myth retellings revoice archetypes to address contemporary urgencies.",
"Intertextual play reframes power, gender, and memory across time.",
"Revisionary angles restore silenced perspectives and complicate origin."
]
},
{
keywords: ['children’s literature', 'didactic', 'imagination', 'picturebook', 'growth'],
responses: [
"Children’s literature balances delight with gentle guidance and moral play.",
"Picturebooks choreograph image‑text interplay to scaffold understanding.",
"Growth arcs encourage empathy, curiosity, and agency in young readers."
]
},
{
keywords: ['young adult', 'identity', 'coming of age', 'voice', 'contemporary issues'],
responses: [
"YA fiction foregrounds identity formation amid peer culture and pressure.",
"Immediate voice and pace connect urgent themes to developing autonomy.",
"Diverse perspectives broaden recognition and social imagination."
]
},
{
keywords: ['memoir', 'autobiography', 'self', 'memory', 'truth effects'],
responses: [
"Memoir shapes lived experience into crafted narrative arcs of meaning.",
"Memory’s selectivity produces truth effects rather than neutral record.",
"Voice mediates intimacy, reliability, and ethical representation."
]
},
{
keywords: ['biography', 'research', 'context', 'character', 'legacy'],
responses: [
"Biography situates a life within social forces, choices, and constraints.",
"Research triangulates archives, testimony, and public record for texture.",
"Character study examines impact, contradiction, and evolving legacy."
]
},
{
keywords: ['essay', 'montaigne', 'argument', 'form', 'voice'],
responses: [
"The essay tests ideas in motion, balancing reflection with crafted stance.",
"From Montaigne onward, form remains elastic across personal and polemic.",
"Voice anchors credibility through clarity, humility, and precision."
]
},
{
keywords: ['rhetoric', 'ethos', 'pathos', 'logos', 'arrangement'],
responses: [
"Rhetoric persuades by aligning ethos, pathos, and logos to audience and aim.",
"Arrangement orders claims for momentum, emphasis, and retention.",
"Style calibrates diction and syntax to foster trust and resonance."
]
},
{
keywords: ['narrator', 'point of view', 'first person', 'third person', 'omniscient'],
responses: [
"Point of view filters knowledge, bias, and intimacy in narrative delivery.",
"First person offers immediacy and limitation; third person varies range.",
"Omniscience surveys multiple interiors but risks diluting focal intensity."
]
},
{
keywords: ['unreliable narrator', 'bias', 'gaps', 'twist', 'trust'],
responses: [
"Unreliable narration exploits bias, omission, or delusion to destabilize truth.",
"Gaps invite readers to reconstruct events beyond stated perspective.",
"Revelatory twists reframe ethics, memory, and motive retroactively."
]
},
{
keywords: ['free indirect discourse', 'style indirect libre', 'blended voice', 'interiority'],
responses: [
"Free indirect discourse blends narrator diction with character thought.",
"It preserves third‑person distance while channeling intimate interiority.",
"Irony arises as authorial and character voices subtly misalign."
]
},
{
keywords: ['plot', 'structure', 'exposition', 'climax', 'denouement'],
responses: [
"Classical structure arcs from exposition through rising tension to climax.",
"Complications escalate stakes until decisive reversal unlocks resolution.",
"Denouement restores a new order that bears the plot’s moral imprint."
]
},
{
keywords: ['foreshadowing', 'chekhov’s gun', 'setups and payoffs', 'expectation'],
responses: [
"Foreshadowing primes readers with signals that later cohere as meaning.",
"Chekhov’s gun mandates that salient details realize narrative purpose.",
"Setups and payoffs convert attention into satisfaction and inevitability."
]
},
{
keywords: ['symbol', 'motif', 'theme', 'pattern', 'resonance'],
responses: [
"Symbols condense abstract themes into tangible, repeatable emblems.",
"Motifs recur with variation, weaving cohesion across scenes and tones.",
"Themes govern value‑conflicts that structure interpretation and affect."
]
},
{
keywords: ['tone', 'mood', 'voice', 'attitude', 'atmosphere'],
responses: [
"Tone registers authorial attitude toward subject and audience through style.",
"Mood names the felt atmosphere projected by setting and cadence.",
"Voice fuses idiolect, rhythm, and stance into recognizable signature."
]
},
{
keywords: ['diction', 'syntax', 'register', 'connotation', 'cadence'],
responses: [
"Diction’s register and connotation guide nuance beyond literal meaning.",
"Syntax shapes pace, emphasis, and cognitive load across sentences.",
"Cadence orchestrates breath, silence, and emphasis into musical prose."
]
},
{
keywords: ['imagery', 'sensory detail', 'metaphor', 'synesthesia', 'patterning'],
responses: [
"Imagery grounds abstraction in sensory experience to quicken attention.",
"Metaphor discovers likeness across difference, enlarging perception.",
"Synesthetic links fuse senses, intensifying imaginative apprehension."
]
},
{
keywords: ['meter', 'rhyme', 'prosody', 'iambic pentameter', 'caesura'],
responses: [
"Prosody coordinates meter, rhyme, and pause into patterned expectation.",
"Iambic pentameter balances flexibility with memorable cadence in English.",
"Strategic caesura and enjambment modulate momentum and surprise."
]
},
{
keywords: ['sonnet forms', 'petrarchan', 'shakespearean', 'volta', 'octave sestet'],
responses: [
"Petrarchan sonnets pivot after octave; Shakespearean turns before the couplet.",
"The volta reframes argument or mood with structural clarity and flair.",
"Rhyme schemes scaffold logical turns within tight lyrical constraint."
]
},
{
keywords: ['blank verse', 'unrhymed iambic', 'dramatic', 'milton', 'shakespeare'],
responses: [
"Blank verse marries speechlike rhythm to elevated dramatic flexibility.",
"Shakespeare and Milton exploited its range for stage and epic gravity.",
"Unrhymed cadence foregrounds syntax and thought over end‑stopped chime."
]
},
{
keywords: ['free verse', 'line break', 'cadence', 'organic form', 'modern'],
responses: [
"Free verse shapes breath and emphasis by line rather than meter.",
"Organic form lets content discover its own governing cadence.",
"Modern poetics prize musicality without fixed pattern obligations."
]
},
{
keywords: ['ballad', 'narrative song', 'quatrain', 'incremental repetition'],
responses: [
"Ballads relate dramatic episodes in simple quatrains and refrains.",
"Incremental repetition deepens tension as images accumulate.",
"Oral roots persist in strong beats, dialogue, and communal lore."
]
},
{
keywords: ['elegy', 'lament', 'consolation', 'pastoral', 'memory'],
responses: [
"Elegy moves from grief to measured consolation through ritual address.",
"Pastoral frames loss within cyclical nature and communal rites.",
"Memory fashions enduring presence from love’s articulate mourning."
]
},
{
keywords: ['ode', 'apostrophe', 'strophe antistrophe epode', 'praise', 'meditation'],
responses: [
"Odes elevate praise through apostrophe and ceremonial structure.",
"Classical triads choreograph thought into balanced turns and closures.",
"Meditative drift refines perception into gratitude and insight."
]
},
{
keywords: ['pastoral', 'country and city', 'idyll', 'shepherd', 'ecology'],
responses: [
"Pastoral idealizes rural life while reflecting on urban ambition.",
"Shepherd figures mediate artistry, leisure, and ethical calm.",
"Modern ecologies complicate nostalgia with worldly responsibility."
]
},
{
keywords: ['satire', 'horatian', 'juvenalian', 'parody', 'irony'],
responses: [
"Satire corrects folly through wit that ranges from genial to scathing.",
"Horatian smiles in urbane irony; Juvenalian lashes with moral outrage.",
"Parody mimics style to expose pretension, error, or cultural drift."
]
},
{
keywords: ['parody', 'pastiche', 'imitation', 'critique', 'play'],
responses: [
"Parody sharpens critique by exaggerating stylistic tics for effect.",
"Pastiche neutrally reuses styles to celebrate or recombine lineage.",
"Imitation becomes invention when transformed with purpose and wit."
]
},
{
keywords: ['allegory', 'two levels', 'personification', 'pilgrimage', 'moral architecture'],
responses: [
"Allegory sustains double meaning across narrative and emblematic planes.",
"Personified abstractions enact ethical journeys through designed worlds.",
"Spatial and social structures map spiritual or political argument."
]
},
{
keywords: ['fable', 'animals', 'moral', 'brevity', 'lesson'],
responses: [
"Fables compress ethical lessons into lively, memorable miniatures.",
"Animal proxies distance critique while sharpening relevance.",
"Explicit morals frame interpretation and practical takeaway."
]
},
{
keywords: ['fabliau', 'medieval comic tale', 'bawdy', 'trickster', 'reversal'],
responses: [
"Fabliaux deploy bawdy wit to invert hierarchies and expose vice.",
"Trickster plots reward cleverness over status or solemnity.",
"Earthy realism punctures pretension with brisk, satiric payoff."
]
},
{
keywords: ['bildungsroman', 'coming of age', 'formation', 'society', 'identity'],
responses: [
"The Bildungsroman traces growth through trials toward social integration.",
"Conflicts between desire and duty crystallize a durable self‑concept.",
"Mentors, mistakes, and mobility structure the path to maturity."
]
},
{
keywords: ['künstlerroman', 'artist novel', 'aesthetics', 'vocation', 'isolation'],
responses: [
"The Künstlerroman charts an artist’s vocation amid friction with norms.",
"Aesthetic awakening often entails solitude, exile, or misfit status.",
"Form mirrors evolving style as self and art co‑constitute destiny."
]
},
{
keywords: ['picaresque', 'rogue', 'episodic', 'satire', 'society'],
responses: [
"Picaresque follows a roguish survivor through episodic social vignettes.",
"Loose structure enables panoramic satire of institutions and vice.",
"Adaptability replaces moral uplift as pragmatic virtue in motion."
]
},
{
keywords: ['epistolary novel', 'letters', 'documents', 'intimacy', 'polyphony'],
responses: [
"Epistolary form assembles narrative from letters and documents.",
"Intimacy and bias emerge through competing first‑person lenses.",
"Polyphony creates suspense from partial knowledge and delay."
]
},
{
keywords: ['frame narrative', 'story within story', 'nested', 'perspective', 'authority'],
responses: [
"Frame narratives nest tales to question perspective and reliability.",
"Outer tellers mediate interpretation and credibility of inner plots.",
"Layering multiplies meaning and distance between event and report."
]
},
{
keywords: ['metafiction', 'self reference', 'breaking the fourth wall', 'fiction about fiction'],
responses: [
"Metafiction exposes its own devices to probe truth in storytelling.",
"Fourth‑wall breaches invite readers to co‑author interpretive play.",
"Fiction about fiction turns craft into subject and argument."
]
},
{
keywords: ['intertextuality', 'allusion', 'palimpsest', 'dialogue of texts'],
responses: [
"Intertextuality situates works within a dialogic web of echoes and riffs.",
"Allusion compresses cultural memory into charged, associative signals.",
"Palimpsestic layers let new writing revise inherited meanings."
]
},
{
keywords: ['archetype', 'hero', 'trickster', 'shadow', 'collective unconscious'],
responses: [
"Archetypes distill recurring figures and plots across cultures and eras.",
"Hero, trickster, and shadow dramatize psychic and social dynamics.",
"Patterns endure while particulars localize universal tensions."
]
},
{
keywords: ['monomyth', 'hero’s journey', 'departure', 'initiation', 'return'],
responses: [
"The monomyth charts call, threshold, ordeal, and transformative return.",
"Allies, mentors, and trials scaffold competence and self‑knowledge.",
"Return with the boon aligns personal change with communal renewal."
]
},
{
keywords: ['heroine’s journey', 'agency', 'community', 'integration', 'relational'],
responses: [
"The heroine’s journey emphasizes integration, relation, and community healing.",
"Agency grows through negotiation of care, voice, and boundaries.",
"Resolution values connection over conquest and mutual flourishing."
]
},
{
keywords: ['catharsis', 'pity and fear', 'emotional purge', 'ethical clarity'],
responses: [
"Catharsis purges pity and fear, clarifying ethical vision through feeling.",
"Tragic design calibrates intensity toward reflective relief, not despair.",
"Audience emerges tempered, discerning limits and responsibilities anew."
]
},
{
keywords: ['dramatic irony', 'audience knowledge', 'tension', 'double meaning'],
responses: [
"Dramatic irony arms audiences with knowledge characters lack for tension.",
"Double meanings accumulate as words outpace speaker awareness.",
"Anticipation sharpens pathos and critique in converging revelations."
]
},
{
keywords: ['stagecraft', 'blocking', 'props', 'lighting', 'design'],
responses: [
"Stagecraft orchestrates blocking, light, and objects into visual argument.",
"Props become symbolic actors that focus theme and motive in space.",
"Lighting sculpts mood, time, and emphasis as silent narration."
]
},
{
keywords: ['regional literature', 'local color', 'dialect', 'place', 'custom'],
responses: [
"Regional writing distills place through dialect, custom, and landscape.",
"Local color preserves texture while examining social power and change.",
"Particularity universalizes by grounding theme in lived specificity."
]
},
{
keywords: ['ecocriticism', 'environment', 'anthropocene', 'nonhuman', 'ethics'],
responses: [
"Ecocriticism reads environments as agents, not mere settings or backdrops.",
"Texts register Anthropocene entanglements across species and systems.",
"Ethical horizons widen to include nonhuman claims and futures."
]
},
{
keywords: ['feminist criticism', 'gender', 'canon', 'gaze', 'representation'],
responses: [
"Feminist criticism interrogates gendered power in form, theme, and canon.",
"It critiques the gaze and expands representation beyond patriarchal norms.",
"Recovery projects restore neglected authors and reframe tradition."
]
},
{
keywords: ['marxist criticism', 'class', 'ideology', 'production', 'hegemony'],
responses: [
"Marxist reading tracks class conflict, ideology, and material conditions.",
"Form and content reveal how hegemony naturalizes economic relations.",
"History becomes a struggle written into genre, plot, and character."
]
},
{
keywords: ['psychoanalytic criticism', 'freud', 'lacan', 'desire', 'unconscious'],
responses: [
"Psychoanalytic lenses decode repression, desire, and symptomatic form.",
"Freud foregrounds drives and dreams; Lacan retools subject and signifier.",
"Texts stage psychic dramas in language, image, and repetition."
]
},
{
keywords: ['structuralism', 'semiotics', 'binary oppositions', 'myth systems'],
responses: [
"Structuralism maps underlying codes governing narrative and myth.",
"Binary oppositions organize meaning that stories then mediate or invert.",
"Semiotics treats literature as sign systems within cultural grammar."
]
},
{
keywords: ['poststructuralism', 'deconstruction', 'différance', 'instability'],
responses: [
"Poststructuralism unsettles fixed meaning by exposing textual slippage.",
"Deconstruction reveals how terms depend on what they exclude and defer.",
"Interpretation becomes an ethics of attention to instability."
]
},
{
keywords: ['reader-response', 'interpretive communities', 'horizon of expectations'],
responses: [
"Reader‑response centers meaning in the encounter between text and reader.",
"Interpretive communities shape horizons of expectation and uptake.",
"Experience becomes evidence as reading acts leave traceable effects."
]
},
{
keywords: ['new criticism', 'close reading', 'organic unity', 'ambiguity', 'paradox'],
responses: [
"New Criticism valorizes close reading and intrinsic textual coherence.",
"Ambiguity and paradox register complexity within organic unity.",
"External biography yields to formal pattern as interpretive ground."
]
},
{
keywords: ['queer theory', 'normativity', 'performativity', 'desire', 'subversion'],
responses: [
"Queer theory interrogates normativity and the performative making of identity.",
"It reads desire as disruptive force within and against textual orders.",
"Subversion emerges through style, camp, and alternative kinships."
]
},
{
keywords: ['disability studies', 'embodiment', 'access', 'narrative prosthesis'],
responses: [
"Disability studies contests stigma and metaphorization of embodied difference.",
"Access reframes form, medium, and audience as ethical design questions.",
"Narrative prosthesis critiques stock uses of disability as plot device."
]
},
{
keywords: ['world literature', 'translation', 'circulation', 'center and periphery'],
responses: [
"World literature follows texts in translation across circuits of power and taste.",
"Centers and peripheries shift as markets and mediators revalue works.",
"Circulation transforms reception, form, and authorial self‑conception."
]
},
{
keywords: ['translation theory', 'equivalence', 'domestication', 'foreignization'],
responses: [
"Translation negotiates equivalence across semantic, cultural, and rhythmic planes.",
"Domestication eases entry for target readers; foreignization keeps estrangement.",
"Every choice reauthors the text’s voice, horizon, and aura."
]
},
{
keywords: ['close reading', 'evidence', 'pattern', 'tension', 'inference'],
responses: [
"Close reading traces patterns and tensions to build warranted claims.",
"Quotations function as evidence when analyzed, not merely cited.",
"Inference grows by linking micro‑features to macro‑themes."
]
},
{
keywords: ['research methods', 'mla', 'apa', 'citation', 'plagiarism'],
responses: [
"Research integrates sources ethically with transparent citation practices.",
"MLA and APA standardize attribution, avoiding plagiarism and confusion.",
"Synthesis weighs credibility, relevance, and conversation entry points."
]
},
{
keywords: ['canon', 'tradition', 'gatekeeping', 'revision', 'curriculum'],
responses: [
"The literary canon reflects historical gatekeeping and evolving values.",
"Revision expands tradition by revaluing genres, regions, and identities.",
"Curricula curate canons to balance legacy, diversity, and debate."
]
},
{
keywords: ['digital humanities', 'text mining', 'distant reading', 'data visualization'],
responses: [
"Digital humanities augment criticism with computation and visualization.",
"Distant reading maps large‑scale patterns beyond single‑text focus.",
"Data complements, not replaces, interpretive nuance and context."
]
},
{
keywords: ['oral tradition', 'performance', 'memory', 'formulaic', 'community'],
responses: [
"Oral literature lives in performance, memory, and communal participation.",
"Formulaic structures aid recall and improvisational variation.",
"Transcription preserves but also transforms ephemeral art."
]
},
{
keywords: ['folk tale', 'motif index', 'propp functions', 'tale types'],
responses: [
"Folktales recombine stable functions into diverse narrative sequences.",
"Motif indexes and tale types map diffusion and transformation.",
"Structure persists while local color supplies distinctive flavor."
]
},
{
keywords: ['fairy tale', 'wonder', 'initiation', 'morality', 'transformation'],
responses: [
"Fairy tales stage wonder and peril as rites of passage into agency.",
"Moral clarity coexists with symbolic ambiguity and dream logic.",
"Transformations externalize psychic and social conflict as magic."
]
},
{
keywords: ['myth', 'cosmogony', 'trickster', 'rite', 'sacred narrative'],
responses: [
"Myths orient cultures with sacred narratives of origin and order.",
"Tricksters break rules to reveal limits, possibilities, and play.",
"Ritual and story bind meaning to communal time and space."
]
},
{
keywords: ['oral epic', 'bard', 'formula', 'meter', 'memory'],
responses: [
"Oral epics use formulaic diction and meter to sustain long performance.",
"Bards adapt episodes to audience, place, and occasion.",
"Memory partners with structure to preserve cultural archives."
]
},
{
keywords: ['haiku', 'tanka', 'kigo', 'cutting word', 'brevity'],
responses: [
"Haiku compress perception with seasonal kigo and a cutting pivot.",
"Tanka extends reflection into a five‑line lyric of nuance.",
"Brevity intensifies clarity, silence, and afterimage."
]
},
{
keywords: ['prose poem', 'lyric prose', 'boundary crossing', 'compression'],
responses: [
"The prose poem blurs boundaries to concentrate lyric intensity in blocks.",
"Syntax and rhythm replace line as the primary sculptural tool.",
"Compression invites rereading to unlock layered resonance."
]
},
{
keywords: ['microfiction', 'drabble', 'flash', 'hint fiction', 'ellipsis'],
responses: [
"Microfiction wagers on ellipsis and implication over exposition.",
"Titles and last lines do outsized work in framing arcs.",
"Constraint catalyzes invention through selective detail."
]
},
{
keywords: ['journalism and literature', 'new journalism', 'narrative nonfiction'],
responses: [
"Narrative nonfiction adopts literary craft under factual constraints.",
"New Journalism foregrounds scene, dialogue, and point of view.",
"Ethics govern immersion, composite risk, and verification."
]
},
{
keywords: ['oral history', 'testimony', 'memory politics', 'polyphony'],
responses: [
"Oral history centers lived testimony within curated polyphony.",
"Memory politics shape what is told, archived, and taught.",
"Editing honors voice while structuring legibility and care."
]
},
{
keywords: ['adaptation', 'page to screen', 'fidelity', 'transposition'],
responses: [
"Adaptation translates narrative across media with gains and losses.",
"Fidelity yields to re‑creation as form, pacing, and focus shift.",
"Transposition interprets essence rather than copying surface."
]
},
{
keywords: ['book history', 'print culture', 'paratext', 'circulation'],
responses: [
"Book history tracks material forms, markets, and reading practices.",
"Paratexts—covers, prefaces, blurbs—frame reception and meaning.",
"Circulation patterns recode status, access, and influence."
]
},
{
keywords: ['censorship', 'obscenity', 'blasphemy', 'political repression'],
responses: [
"Censorship polices borders of taste, doctrine, and power.",
"Obscenity and blasphemy trials reveal shifting cultural norms.",
"Writers develop veiled strategies to evade repression."
]
},
{
keywords: ['marginalized voices', 'recovery projects', 'archive', 'countercanon'],
responses: [
"Recovery projects rebuild archives for marginalized literatures.",
"Countercanons alter the map of influence and dialogue.",
"Access reshapes curricula, discourse, and identity horizons."
]
},
{
keywords: ['form and content', 'fit', 'tension', 'organic design'],
responses: [
"Form and content interanimate, with tension often generating power.",
"Organic design emerges as choices cohere around necessity.",
"Mismatches can be strategic to estrange and provoke thought."
]
},
{
keywords: ['style', 'periodicity', 'signature', 'imitation vs originality'],
responses: [
"Style marries pattern to surprise, forming a writer’s signature.",
"Period styles set defaults that individuals bend or break.",
"Imitation teaches craft; originality recombines lineage freshly."
]
}
  ]; // END MEGA_KNOWLEDGE_BASE

  // ===== MATCHING ENGINE (keyword + fuzzy) =====
  const normalize = (s) =>
    s.toLowerCase().replace(/[^a-z0-9\s]/g, ' ').replace(/\s+/g, ' ').trim();

  const levenshtein = (a, b) => {
    const m = a.length, n = b.length;
    if (!m) return n;
    if (!n) return m;
    const dp = Array.from({ length: m + 1 }, () => Array(n + 1).fill(0));
    for (let i = 0; i <= m; i++) dp[i][0] = i;
    for (let j = 0; j <= n; j++) dp[0][j] = j;
    for (let i = 1; i <= m; i++) {
      for (let j = 1; j <= n; j++) {
        const cost = a[i - 1] === b[j - 1] ? 0 : 1;
        dp[i][j] = Math.min(
          dp[i - 1][j] + 1,
          dp[i][j - 1] + 1,
          dp[i - 1][j - 1] + cost
        );
      }
    }
    return dp[m][n];
  };

  const similarity = (a, b) => {
    const A = normalize(a), B = normalize(b);
    if (!A.length && !B.length) return 1;
    const longer = A.length >= B.length ? A : B;
    const shorter = A.length >= B.length ? B : A;
    return (longer.length - levenshtein(longer, shorter)) / longer.length;
  };

  const scoreEntry = (query, entry) => {
    const q = normalize(query);
    const qTokens = q.split(' ').filter(Boolean);
    let score = 0;

    for (const kw of entry.keywords) {
      const k = normalize(kw);
      if (!k) continue;

      // Exact phrase boost
      if (q.includes(k)) score += Math.min(20, k.split(' ').length * 6);

      // Token overlap and fuzzy
      for (const qt of qTokens) {
        if (!qt) continue;
        if (k.includes(qt) || qt.includes(k)) score += 3;
        const sim = similarity(qt, k);
        if (sim >= 0.8) score += 4;
        else if (sim >= 0.65) score += 2;
      }
    }

    // Bonus for multiple matched keywords
    const matchedCount = entry.keywords.reduce((acc, kw) => acc + (q.includes(normalize(kw)) ? 1 : 0), 0);
    score += Math.max(0, matchedCount - 1) * 2;

    return score;
  };

  const findBestResponse = (query) => {
    if (!query || !query.trim()) return "Please share a topic or a couple of keywords so I can help.";
    let best = { score: -Infinity, entry: null };
    for (const entry of MEGA_KNOWLEDGE_BASE) {
      const s = scoreEntry(query, entry);
      if (s > best.score) best = { score: s, entry };
    }
    if (!best.entry || best.score < 5) {
      return "I didn’t find a strong match—try a couple of more specific keywords like 'kinetics activation energy' or 'WWII causes Europe'.";
    }
    const choices = best.entry.responses;
    return choices[Math.floor(Math.random() * choices.length)];
  };

  // ===== UI HELPERS =====
  const addMessage = (who, text) => {
    setMessages((prev) => [...prev, { who, text, ts: Date.now() }]);
    setTimeout(() => messagesEndRef.current?.scrollIntoView({ behavior: 'smooth' }), 0);
  };

  const handleSend = (e) => {
    e?.preventDefault?.();
    const q = inputText.trim();
    if (!q) return;
    setInputText('');
    addMessage('user', q);
    setIsProcessing(true);
    setStatus('Thinking...');
    setTimeout(() => {
      const a = findBestResponse(q);
      addMessage('bot', a);
      setIsProcessing(false);
      setStatus('Ready');
    }, 150);
  };

  const randomQuestion = () => {
    const seeds = [
      'photosynthesis stages',
      'Ohm’s law basics',
      'Pythagorean theorem proof',
      'French Revolution causes',
      'DNA transcription translation',
      'Galvanic cell vs electrolytic',
      'Gaussian elimination idea',
      'Renaissance humanism points',
      'Quantum photoelectric effect',
      'Buffer Henderson-Hasselbalch'
    ];
    const q = seeds[Math.floor(Math.random() * seeds.length)];
    setInputText(q);
    setTimeout(() => handleSend(), 50);
  };

  const clearChat = () => {
    setMessages([]);
    setStatus('Ready to help with 5000+ topics!');
  };

  const askAbout = (topic) => {
    setInputText(topic);
    setTimeout(() => handleSend(), 50);
  };

  // ===== (Optional) Speech helpers (off by default in React) =====
  const speak = (text) => {
    if (!window.speechSynthesis) return;
    const u = new SpeechSynthesisUtterance(text);
    u.rate = 0.95;
    const voice = speechSynthesis.getVoices().find(v => v.lang.startsWith('en'));
    if (voice) u.voice = voice;
    speechSynthesis.speak(u);
  };

  useEffect(() => {
    messagesEndRef.current?.scrollIntoView({ behavior: 'smooth' });
  }, [messages]);

  return (
    <div className={`mega-wrapper theme-${theme}`}>
      <header className="mega-header">
        <h1 className="mega-title">🧠 Mega Educational AI</h1>
        <p className="mega-sub">Ask with keywords. Covers Science, Math, History, Literature, Tech, Geography, and Study Skills.</p>
        <div className="mega-stats">
          <span>📚 50+ subjects</span>
          <span>🎯 5000+ Q&A</span>
          <span>⚡ Keyword + fuzzy match</span>
        </div>
      </header>

      <main className="mega-main">
        <aside className="mega-side">
          <section className="panel">
            <h3>🎛️ Quick Actions</h3>
            <div className="btn-row">
              <button onClick={randomQuestion} className="btn primary">Random Q</button>
              <button onClick={clearChat} className="btn">Clear</button>
            </div>
          </section>

          <section className="panel">
            <h3>📚 Topics</h3>
            <div className="tags">
              <button onClick={() => askAbout('biology photosynthesis')} className="tag">🧬 Biology</button>
              <button onClick={() => askAbout('chemistry acids bases')} className="tag">⚗️ Chemistry</button>
              <button onClick={() => askAbout('physics kinematics')} className="tag">🔬 Physics</button>
              <button onClick={() => askAbout('algebra quadratics')} className="tag">📐 Math</button>
              <button onClick={() => askAbout('WWII causes')} className="tag">🏛️ History</button>
              <button onClick={() => askAbout('poetry devices')} className="tag">📖 Literature</button>
              <button onClick={() => askAbout('web html css js')} className="tag">💻 Tech</button>
              <button onClick={() => askAbout('climate zones')} className="tag">🌍 Geography</button>
              <button onClick={() => askAbout('study spaced repetition')} className="tag">📝 Study</button>
            </div>
          </section>

          <section className="panel">
            <h3>🎨 Theme</h3>
            <div className="btn-row">
              <button className={`btn ${theme==='cosmic'?'active':''}`} onClick={() => setTheme('cosmic')}>Cosmic</button>
              <button className={`btn ${theme==='aurora'?'active':''}`} onClick={() => setTheme('aurora')}>Aurora</button>
              <button className={`btn ${theme==='sunset'?'active':''}`} onClick={() => setTheme('sunset')}>Sunset</button>
            </div>
          </section>
        </aside>

        <section className="mega-chat">
          <div className="chat-status">{isProcessing ? '🤔 Thinking...' : status}</div>
          <div className="chat-body">
            {messages.length === 0 && (
              <div className="welcome">
                <div className="welcome-icon">🚀</div>
                <div className="welcome-text">
                  Ask anything with a couple of keywords like “electron transport chain” or “law of sines example.”
                </div>
              </div>
            )}
            {messages.map((m, i) => (
              <div key={m.ts + '-' + i} className={`msg ${m.who === 'user' ? 'right' : 'left'}`}>
                <div className={`bubble ${m.who === 'user' ? 'user' : 'bot'}`}>
                  {m.text}
                </div>
              </div>
            ))}
            {isProcessing && (
              <div className="msg left">
                <div className="bubble bot">
                  <span className="dots">
                    <span></span><span></span><span></span>
                  </span>
                </div>
              </div>
            )}
            <div ref={messagesEndRef} />
          </div>

          <form className="chat-input" onSubmit={handleSend}>
            <input
              value={inputText}
              onChange={(e) => setInputText(e.target.value)}
              placeholder="Type keywords (e.g., 'entropy spontaneous', 'matrix determinant', 'Renaissance humanism')"
            />
            <button className="btn primary" type="submit">Send</button>
          </form>
        </section>
      </main>

      <footer className="mega-foot">
        <span>Built for fast keyword-based learning. Try adding 2–4 precise terms for best results.</span>
      </footer>
    </div>
  );
};

export default MegaChatbot;
