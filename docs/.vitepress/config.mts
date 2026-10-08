import { withMermaid } from 'vitepress-plugin-mermaid'
import apiSidebar from '../api/_sidebar.json'

// one sidebar for every page but the API, in reading order: the "previous" and "next" links at the bottom of a page
// follow it, so the documentation reads from the first page to the last
const docs = [
  {
    text: 'Start here',
    items: [
      { text: 'Introduction', link: '/guide/' },
      { text: 'Installation', link: '/guide/installation' },
      { text: 'Quick start',  link: '/guide/quickstart' },
    ],
  },
  {
    text: 'Schemes',
    items: [
      { text: 'Catalogue',           link: '/guide/schemes' },
      { text: 'Runge-Kutta schemes', link: '/guide/schemes-runge-kutta' },
      { text: 'Multistep schemes',   link: '/guide/schemes-multistep' },
    ],
  },
  {
    text: 'Verification',
    items: [
      { text: 'Testing', link: '/guide/testing' },
    ],
  },
  {
    text: 'Examples',
    items: [
      { text: 'Oscillation equations', link: '/guide/example-oscillation' },
      { text: '1D Euler equations',    link: '/guide/example-euler-1d' },
    ],
  },
  {
    text: 'Project',
    items: [
      { text: 'Changelog',    link: '/guide/changelog' },
      { text: 'Contributing', link: '/guide/contributing' },
    ],
  },
]

export default withMermaid({
  title: 'FOODIE',
  description: 'Fortran Object-Oriented Differential-equations Integration Environment',
  base: '/FOODIE/',

  markdown: {
    math: true,
    languages: ['fortran-free-form', 'fortran-fixed-form'],
    languageAlias: {
      fortran: 'fortran-free-form',
      f90: 'fortran-free-form',
      f03: 'fortran-free-form',
      f08: 'fortran-free-form',
    },
  },

  themeConfig: {
    nav: [
      { text: 'Home', link: '/' },
      { text: 'Start here', link: '/guide/', activeMatch: '^/guide/(index|installation|quickstart)' },
      { text: 'Schemes', link: '/guide/schemes', activeMatch: '^/guide/schemes' },
      { text: 'Testing', link: '/guide/testing', activeMatch: '^/guide/testing' },
      { text: 'Examples', link: '/guide/example-oscillation', activeMatch: '^/guide/example-' },
      { text: 'API', link: '/api/' },
      {
        text: 'Project',
        items: [
          { text: 'Changelog',    link: '/guide/changelog' },
          { text: 'Contributing', link: '/guide/contributing' },
          { text: 'GitHub',       link: 'https://github.com/Fortran-FOSS-Programmers/FOODIE' },
        ],
      },
    ],

    sidebar: {
      '/guide/': docs,
      '/api/': [
        {
          text: 'API Reference',
          items: [
            { text: 'Overview', link: '/api/' },
          ],
        },
        ...apiSidebar,
      ],
    },

    socialLinks: [
      { icon: 'github', link: 'https://github.com/Fortran-FOSS-Programmers/FOODIE' },
    ],

    search: {
      provider: 'local',
    },

    footer: {
      message: 'Released under the GPL v3, BSD 2-Clause, BSD 3-Clause and MIT licenses.',
      copyright: 'Copyright © 2026 Stefano Zaghi',
    },
  },

  mermaid: {},

  vite: {
    // Build with an explicit modern JS target so the docs compile regardless of
    // which mermaid/vitepress/esbuild versions npm resolves. Vite's default
    // es2020 target forces esbuild to down-level modern syntax (e.g. the
    // destructuring mermaid 11.16+ emits), which it refuses to do and the build
    // dies. es2022 needs no lowering and is within VitePress's browser floor.
    build: {
      target: 'es2022',
    },
  },
})
