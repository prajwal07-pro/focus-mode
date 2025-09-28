import React, { useEffect, useState } from 'react';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip } from 'recharts';
import { subscribeFocus } from '../services/analytics';

export default function Analytics({ sessionId }) {
  const [data, setData] = useState([]);
  useEffect(() => {
    let i = 0;
    const unsub = subscribeFocus(sessionId, (row) => {
      setData((d) => [...d, { timestamp: ++i, focus_score: row.focus_score }]);
    });
    return () => unsub && unsub();
  }, [sessionId]);

  return (
    <LineChart width={500} height={300} data={data}>
      <CartesianGrid stroke="#2a3252" />
      <XAxis dataKey="timestamp" />
      <YAxis />
      <Tooltip />
      <Line type="monotone" dataKey="focus_score" stroke="#8884d8" />
    </LineChart>
  );
}
