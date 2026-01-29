// @ts-check
// `@type` JSDoc annotations allow editor autocompletion and type checking
// (when paired with `@ts-check`).
// There are various equivalent ways to declare your Docusaurus config.
// See: https://docusaurus.io/docs/api/docusaurus-config

import {themes as prismThemes} from 'prism-react-renderer';

// This runs in Node.js - Don't use client-side code here (browser APIs, JSX...)

/** @type {import('@docusaurus/types').Config} */
const config = {
  title: 'TUM Data Knowledge Hub',
  tagline: 'Cultivating Research Integrity and Reproducible Science',
  favicon: 'img/tum-logo.svg',

  // Future flags, see https://docusaurus.io/docs/api/docusaurus-config#future
  future: {
    v4: true, // Improve compatibility with the upcoming Docusaurus v4
  },

  // Set the production url of your site here
  url: 'https://tum-research-data-hub.github.io',
  // Set the /<baseUrl>/ pathname under which your site is served
  // For GitHub pages deployment, it is often '/<projectName>/'
  baseUrl: '/knowledge-hub/',

  // GitHub pages deployment config.
  // If you aren't using GitHub pages, you don't need these.
  organizationName: 'tum-research-data-hub', // Usually your GitHub org/user name.
  projectName: 'knowledge-hub', // Usually your repo name.
  deploymentBranch: 'gh-pages', // Explicitly tells Docusaurus which branch to push the build to
  trailingSlash: false,  // Helps with URL consistency on GitHub Pages

  onBrokenLinks: 'throw', // Using 'throw' ensures that you never accidentally publish a broken website to your users
  onBrokenMarkdownLinks: 'warn', // This specifically checks internal links inside your Markdown files

  // Even if you don't use internationalization, you can use this field to set
  // useful metadata like html lang. For example, if your site is Chinese, you
  // may want to replace "en" with "zh-Hans".

  i18n: {
    defaultLocale: 'en',
    locales: ['en'],
  },

// This is the section to add a local search functionality 
themes: [
    [
      require.resolve("@easyops-cn/docusaurus-search-local"),
      {
        hashed: true,
        language: ["en"],
        highlightSearchTermsOnTargetPage: true,
        explicitSearchResultPath: true,
        indexDocs: true,
        indexBlog: false,
        docsRouteBasePath: "/",
      },
    ],
  ],

  presets: [
    [
      'classic',
      /** @type {import('@docusaurus/preset-classic').Options} */
      ({
        docs: {
          sidebarPath: './sidebars.js',
          routeBasePath: '/', 
          // Please change this to your repo.
          // Remove this to remove the "edit this page" links.
          editUrl:
            'https://github.com/tum-research-data-hub/knowledge-hub/tree/main/',
        },
        blog: false,
        theme: {
          customCss: './src/css/custom.css',
        },
      }),
    ],
  ],

  themeConfig:
    /** @type {import('@docusaurus/preset-classic').ThemeConfig} */
    ({
      // Replace with your project's social card
      image: 'img/tum-rdhub-social-card.webp',
      colorMode: {
        respectPrefersColorScheme: true,
      },
      navbar: {
        title: 'Data Knowledge Hub',
        logo: {
          alt: 'TUM Logo',
          src: 'img/tum-logo.svg',
        },
        items: [
          {
            type: 'docSidebar',
            sidebarId: 'mySidebar',
            position: 'left',
            label: 'Knowledge Hub',
          },
          {
            href: 'https://github.com/tum-research-data-hub/knowledge-hub',
            label: 'GitHub',
            position: 'right',
          },
          {
            href: 'mailto:researchdata@tum.de',
            label: 'Contact',
            position: 'right',
          },
        ],
      },
      footer: {
        style: 'light',
        links: [
          {
            title: 'Tools',
            items: [
              { label: 'TUM DataTagger', href: 'https://www.ub.tum.de/en/datatagger' },
              { label: 'mediaTUM', href: 'https://mediatum.ub.tum.de/' },
              { label: 'eLabFTW', href: 'https://www.ub.tum.de/en/eln' },
            ],
          },
          {
            title: 'Support',
            items: [
              { label: 'MDSI', href: 'https://www.mdsi.tum.de/' },
              { label: 'TUM University Library', href: 'https://www.ub.tum.de/' },
            ],
          },
        ],
        copyright: `Copyright © ${new Date().getFullYear()} Technical University of Munich. Built with Docusaurus.`,
      },
      prism: {
        theme: prismThemes.github,
        darkTheme: prismThemes.dracula,
      },
    tableOfContents: {
        minHeadingLevel: 2,
        maxHeadingLevel: 4,
      }
    }),
};

export default config;
