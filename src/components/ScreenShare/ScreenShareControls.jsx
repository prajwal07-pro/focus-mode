import React from 'react';
//import styles from './ScreenShareControls.module.css';

const ScreenShareControls = ({ isActive, onStop, onToggleVisibility }) => {
  if (!isActive) return null;

  return (
    <div id="screen-share-controls" className={styles.controls}>
      <div className={styles.indicator}>
        <span className={styles.recordingDot}>🔴</span>
        <span>Sharing screen</span>
      </div>
      <div className={styles.buttons}>
        <button onClick={onToggleVisibility}>Hide</button>
        <button onClick={onStop} className={styles.stopButton}>
          Stop Sharing
        </button>
      </div>
    </div>
  );
};

export default ScreenShareControls;
