import React, { useState } from 'react';
import FocusMode from './components/FocusMode.jsx';
import EducationalChatbot from './components/EducationalChatbot.jsx';

export default function App() {
  const [sessionId, setSessionId] = useState('demo-session');
  
  return (
    <div className="app">
      <main className="adaptive-grid">
        <section className="focus-section-main">
          <h2>🎯 Focus Mode</h2>
          <FocusMode sessionId={sessionId} />
        </section>

        <section className="chatbot-section-side">
          <h2>🤖 Study Assistant</h2>
          <EducationalChatbot />
        </section>
      </main>
    </div>
  );
}
