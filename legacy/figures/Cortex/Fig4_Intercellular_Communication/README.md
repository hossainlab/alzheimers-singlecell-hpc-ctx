# Figure 4: Intercellular Communication Networks (Prefrontal Cortex)

> **Plain-Language Summary for Collaborators**  
> Brain cells do not work in isolation. Astrocytes talk to blood vessels, microglia listen for neuronal damage, and blood vessels receive instructions from pericytes.  
> In Step 4, we use **CellChat v2.2** to eavesdrop on the cellular telephone calls between cell types. We measure which communication lines are working normally, which went dead, and which became toxic and hyper-active in Alzheimer's Disease.

---

## 🖼️ Panel Walkthrough (What You Are Looking At)

### 1. `Fig4a_intercellular_interaction_strength.*` — Total Telephone Call Volume
- **What it shows**: Side-by-side bar plots comparing the total number of cellular interactions (left) and the overall communication strength (right) between Control and Alzheimer's.
- **In plain English**: Shows how the disease reshapes the overall volume of cellular messages passing through the cortical neurovascular unit.

### 2. `Fig4b_differential_interaction_network.*` — Broken vs. Jammed Phone Lines
- **What it shows**: A circular diagram of all cell types showing net changes in communication.
- **How to read the colors**:
  - **Red lines**: Communication channels that are abnormally **hyper-activated** or overloaded in Alzheimer's.
  - **Blue lines**: Healthy communication channels that were **silenced or lost** in Alzheimer's.
  - Line thickness indicates how much the signal changed.

### 3. `Fig4c_signaling_information_flow_ranking.*` — Conversation Topic Ranking
- **What it shows**: A ranked bar chart of 27 molecular communication pathways (topics of conversation).
- **In plain English**: Tells us which conversation topics changed the most. Pathways like **APP (amyloid precursor protein signaling)** and inflammatory pathways surge to the top in Alzheimer's.

### 4. `Fig4d_nvu_pathway_circos_network.*` — The Amyloid Broadcast Network (Circos Wheel)
- **What it shows**: A circular chord diagram specifically tracing *APP* (amyloid) signaling across all cell types in Control (left) vs Alzheimer's (right).
- **Key finding**: In Alzheimer's, astrocytes and neurons broadcast massive amounts of APP ligand directly onto microglia and endothelial cells, driving blood–brain barrier damage and inflammatory activation.

### 5. `Fig4e_ligand_receptor_communication_bubble.*` — The Exact Molecular Words Spoken
- **What it shows**: A high-resolution bubble plot pinpointing specific molecular "keys and locks" (ligand–receptor pairs like *APP–CD46* and *APP–SORL1*).
- **How to read it**:
  - Circle size indicates how confident we are the signal is active ($p$-value).
  - Red color intensity indicates the communication strength/probability.

---

## 📊 Quantified Data Tables in this Folder
- `Cortex_CellChat_interactions_Control.csv`: Complete catalog of every validated ligand–receptor communication link in healthy control cortex.
- `Cortex_CellChat_interactions_AD.csv`: Complete catalog of every validated communication link in Alzheimer's cortex.
