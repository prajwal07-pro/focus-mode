import React, { useEffect, useState } from 'react';
import * as qna from '@tensorflow-models/qna';

export default function QA() {
  const [model, setModel] = useState(null);
  const [context, setContext] = useState('Photosynthesis is the process by which plants convert light energy into chemical energy.');
  const [question, setQuestion] = useState('What is photosynthesis?');
  const [answers, setAnswers] = useState([]);

  useEffect(() => { qna.load().then(setModel); }, []);

  async function ask() {
    if (!model) return;
    const res = await model.findAnswers(question, context);
    setAnswers(res);
  }

  return (
    <div>
      <textarea rows={4} value={context} onChange={(e) => setContext(e.target.value)} style={{ width: '100%' }} />
      <input value={question} onChange={(e) => setQuestion(e.target.value)} style={{ width: '100%', marginTop: 8 }} />
      <button onClick={ask} style={{ marginTop: 8 }}>Ask</button>
      <ul>
        {answers.map((a, idx) => (
          <li key={idx} className="small">{a.text} (score {a.score.toFixed(3)})</li>
        ))}
      </ul>
    </div>
  );
}
