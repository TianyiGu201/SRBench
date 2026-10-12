module.exports = {
  content: ['./index.html', './{pages,components,services}/**/*.{ts,tsx}', './App.tsx'],
  theme: {extend: {}},
  plugins: [require('@tailwindcss/typography')],
};
