import React from 'react';
//import styles from './FocusMode.module.css';

const FocusModeOverlay = ({ onExit, isScreenSharing }) => {
  return (
    <div className={styles.overlay}>
      <div className={styles.overlayContent}>
        <h1>Focus Mode Active</h1>
        {isScreenSharing && (
          <p className={styles.shareIndicator}>
            🔴 Screen sharing continues in background
          </p>
        )}
        <p>Press Ctrl/Cmd + Shift + F to exit</p>
        <button onClick={onExit} className={styles.exitButton}>
          Exit Focus Mode
        </button>
      </div>
    </div>
  );
};

export default FocusModeOverlay;
