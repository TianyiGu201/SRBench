import React from 'react';
import { NavLink, Outlet } from 'react-router-dom';
import { LayoutDashboard, Activity, Grid, Layers, FileText } from 'lucide-react';

const NavItem: React.FC<{ to: string; icon: React.ReactNode; label: string }> = ({ to, icon, label }) => (
  <NavLink
    to={to}
    className={({ isActive }) =>
      `flex items-center px-4 py-3 text-sm font-medium rounded-lg transition-colors ${
        isActive
          ? 'bg-blue-600 text-white'
          : 'text-gray-300 hover:bg-gray-800 hover:text-white'
      }`
    }
  >
    <span className="mr-3">{icon}</span>
    {label}
  </NavLink>
);

const Layout: React.FC = () => {
  return (
    <div className="flex h-screen bg-gray-50 overflow-hidden font-sans">
      {/* Sidebar */}
      <aside className="w-64 bg-gray-900 text-white flex flex-col flex-shrink-0">
        <div className="p-6 border-b border-gray-800">
          <h1 className="text-2xl font-bold tracking-tight bg-gradient-to-r from-blue-400 to-teal-400 bg-clip-text text-transparent">
            SRBench
          </h1>
          <p className="text-xs text-gray-500 mt-1">Exploring Super-Resolution Strategies for Spatial Transcriptomics</p>
        </div>
        
        <nav className="flex-1 px-4 py-6 space-y-2 overflow-y-auto">
          <NavItem to="/study-design" icon={<FileText size={20} />} label="Platform Design" />
          <NavItem to="/dataset-overview" icon={<LayoutDashboard size={20} />} label="Dataset Overview" />
          <NavItem to="/metric-vis" icon={<Activity size={20} />} label="Metric Visualization" />
          <NavItem to="/out-of-sample" icon={<Grid size={20} />} label="Out-of-sample Prediction" />
          <NavItem to="/diff-resolution" icon={<Layers size={20} />} label="Initial Resolution" />         
        </nav>

        <div className="p-4 border-t border-gray-800 text-xs text-gray-500 text-center">
           © 2026 SRBench Project
        </div>
      </aside>

      {/* Main Content */}
      <main className="flex-1 overflow-auto relative">
        <Outlet />
      </main>
    </div>
  );
};

export default Layout;
