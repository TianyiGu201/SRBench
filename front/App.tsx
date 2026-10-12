import React from 'react';
import { HashRouter, Routes, Route, Navigate } from 'react-router-dom';
import Layout from './components/Layout';
import DatasetOverview from './pages/DatasetOverview';
import MetricVisualization from './pages/MetricVisualization';
import OutOfSamplePrediction from './pages/OutOfSamplePrediction';
import DiffResPrediction from './pages/DiffResPrediction';
import StudyDesign from './pages/StudyDesign';

const App: React.FC = () => {
  return (
    <HashRouter>
      <Routes>
        <Route path="/" element={<Layout />}>
          <Route index element={<Navigate to="/study-design" replace />} />
          <Route path="study-design" element={<StudyDesign />} />
          <Route path="dataset-overview" element={<DatasetOverview />} />
          <Route path="metric-vis" element={<MetricVisualization />} />
          <Route path="out-of-sample" element={<OutOfSamplePrediction />} />
          <Route path="diff-resolution" element={<DiffResPrediction />} />
          <Route path="*" element={<Navigate to="/study-design" replace />} />
        </Route>
      </Routes>
    </HashRouter>
  );
};

export default App;
