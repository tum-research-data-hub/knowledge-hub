## 

# **VS Code Project View: TUM Data Knowledge Hub**

## **This file should be saved in the ROOT directory of your project: know-hub-docu/vscode\_structure.md**

## 

## **📂 File Explorer View (Root Directory)**

know-hub-docu/  
├── 📂 docs/                     \# Documentation content folder  
│   ├── 📂 01-about-the-knowledge-hub/  
│   │   ├── \_category\_.json      \# Sets label to "About the Hub"  
│   │   ├── index.md  
│   │   └── organizational.md  
│   ├── 📂 02-general-knowledge/  
│   │   ├── \_category\_.json      \# Sets label to "General Knowledge"  
│   │   ├── data-types.md  
│   │   ├── repositories.md  
│   │   └── ...  
│   ├── 📂 03-domain-specific-knowledge/  
│   │   ├── \_category\_.json      \# Sets label to "Domain Specific"  
│   │   └── 📂 Bioinformatics-Biology/  
│   │       ├── \_category\_.json  \# Sets label to "Bioinformatics & Biology"  
│   │       └── ...  
│   └── 📂 04-use-cases/  
│       ├── \_category\_.json      \# Sets label to "Use Cases"  
│       └── ...  
├── 📂 static/                   \# Static assets (images, logos)  
│   └── 📂 img/  
│       └── tum-logo.svg  
├── 📂 src/                      \# Custom React components and CSS  
│   └── 📂 css/  
│       └── custom.css  
├── docusaurus.config.js         \# Main site configuration  
├── sidebars.js                  \# Sidebar structure  
├── package.json                 \# Project dependencies  
└── vscode\_structure.md          \# \<--- SAVE THIS FILE HERE

## **🛠️ Summary of Benefits**

1. **Safety:** Placing this file in the root ensures Docusaurus ignores it during the build process (it only processes the docs/ folder by default).  
2. **Order:** Folders in docs/ keep their prefixes (01, 02\) so they stay in order in VS Code, while \_category\_.json handles the clean display names.  
3. **Reference:** You can quickly open this file in VS Code to remember the naming conventions or project roadmap.